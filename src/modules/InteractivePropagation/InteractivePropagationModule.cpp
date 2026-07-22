/**
 * @file
 * @brief Implementation of InteractivePropagation module
 *
 * @copyright Copyright (c) 2017-2024 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "InteractivePropagationModule.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <iomanip>
#include <Eigen/Core>

#include "core/utils/distributions.h"
#include "core/utils/log.h"
#include "objects/exceptions.h"
#include "tools/runge_kutta.h"
#include "objects/PropagationSummary.hpp"


using namespace allpix;


//This constructor is copied from TransientPropagation.cpp with a few added lines for a few useful constants
InteractivePropagationModule::InteractivePropagationModule(Configuration& config,
    Messenger* messenger,
    std::shared_ptr<Detector> detector)
: Module(config, detector), messenger_(messenger), detector_(std::move(detector)) {

    // Save detector model
    model_ = detector_->getModel();

    // Require deposits message for single detector:
    messenger_->bindSingle<DepositedChargeMessage>(this, MsgFlags::REQUIRED);

    // Set default value for config variables
    config_.setDefault<double>("timestep", Units::get(0.01, "ns"));
    config_.setDefault<double>("integration_time", Units::get(25, "ns"));
    config_.setDefault<unsigned int>("charge_per_step", 1);
    config_.setDefault<unsigned int>("max_charge_groups", 1000);
    config_.setDefault<double>("coulomb_distance_limit", Units::get(200.0,"um"));
    config_.setDefault<double>("coulomb_field_limit", Units::get(4e5,"V/cm")); 
    // Rickard 2026-04-05: Added bool for outputting propagation summary objects
    config_.setDefault<bool>("output_propagation_summary", false);
    config_.setDefault<double>("output_propagation_summary_step", config_.get<double>("timestep"));

    // Models:
    config_.setDefault<std::string>("mobility_model", "jacoboni");
    config_.setDefault<std::string>("recombination_model", "none");
    config_.setDefault<std::string>("trapping_model", "none");
    config_.setDefault<std::string>("detrapping_model", "none");

    config_.setDefault<double>("temperature", 293.15);
    config_.setDefault<unsigned int>("distance", 1);
    config_.setDefault<bool>("ignore_magnetic_field", false);
    config_.setDefault<double>("relative_permittivity", 1.0);
    config_.setDefault<double>("surface_reflectivity", 0.0);

    // Set defaults for charge carrier multiplication (not used currently)
    config_.setDefault<double>("multiplication_threshold", 1e-2);
    config_.setDefault<unsigned int>("max_multiplication_level", 5);
    config_.setDefault<std::string>("multiplication_model", "none");

    // Set defaults for extra configurability
    config_.setDefault<bool>("enable_diffusion", true);
    config_.setDefault<bool>("enable_coulomb_repulsion", true);

    config_.setDefault<bool>("propagate_electrons", true);
    config_.setDefault<bool>("propagate_holes", true);

    config_.setDefault<bool>("include_mirror_charges", false);

    // Set defaults for plots
    config_.setDefault<bool>("output_linegraphs", false);
    config_.setDefault<bool>("output_linegraphs_collected", false);
    config_.setDefault<bool>("output_linegraphs_recombined", false);
    config_.setDefault<bool>("output_linegraphs_trapped", false);
    config_.setDefault<bool>("output_animations", false);
    config_.setDefault<bool>("output_rms", false);
    config_.setDefault<bool>("output_plots",
    config_.get<bool>("output_linegraphs") || config_.get<bool>("output_animations") || config_.get<bool>("output_rms"));
    config_.setDefault<bool>("output_animations_color_markers", false);
    config_.setDefault<double>("output_plots_step", config_.get<double>("timestep"));
    config_.setDefault<bool>("output_plots_use_pixel_units", false);
    config_.setDefault<bool>("output_plots_align_pixels", false);
    config_.setDefault<double>("output_plots_theta", 0.0f);
    config_.setDefault<double>("output_plots_phi", 0.0f);

    // Rickard 2026-04-05: Get bool for outputting propagation summary objects
    output_propagation_summary_ = config_.get<bool>("output_propagation_summary");

    // Copy some variables from configuration to avoid lookups:
    temperature_ = config_.get<double>("temperature");
    timestep_ = config_.get<double>("timestep");
    integration_time_ = config_.get<double>("integration_time");
    distance_ = config_.get<unsigned int>("distance");
    charge_per_step_ = config_.get<unsigned int>("charge_per_step");
    max_charge_groups_ = config_.get<unsigned int>("max_charge_groups");
    boltzmann_kT_ = Units::get(8.6173333e-5, "eV/K") * temperature_;
    coulomb_K_ =  1.43996454e-12; //Units::get(1.43996454e-12, "MV mm/e");
    surface_reflectivity_ = config_.get<double>("surface_reflectivity");

    max_multiplication_level_ = config.get<unsigned int>("max_multiplication_level");

    enable_diffusion_ = config.get<bool>("enable_diffusion");
    enable_coulomb_repulsion_ = config.get<bool>("enable_coulomb_repulsion");

    propagate_electrons_ = config.get<bool>("propagate_electrons");
    propagate_holes_ = config.get<bool>("propagate_holes");

    include_mirror_charges_ = config.get<bool>("include_mirror_charges");

    relative_permittivity_ = config.get<double>("relative_permittivity"); // the permativity of materials isn't in allpix, so rely on user to pass this in

    if (enable_coulomb_repulsion_ && relative_permittivity_ == 1) {
        LOG(WARNING) << "Coulomb repulsion is enabled but relative permittivity is set to one. Check that the parameter relative_permittivity isn't misspelled or ommited.";
    }

    const auto coulomb_distance_limit = config_.get<double>("coulomb_distance_limit");
    
    coulomb_distance_limit_squared_ = coulomb_distance_limit * coulomb_distance_limit;
    // coulomb_field_limit_ = config.get<double>("coulomb_field_limit"); // Convert from V/cm to MV/mm (internal field units)

    // Store the coulomb_field_limit_ 

    const auto configured_coulomb_field_limit =
    config_.get<double>("coulomb_field_limit");

    coulomb_field_limit_ =
        configured_coulomb_field_limit;

    LOG(WARNING)
        << "[COULOMB_DEBUG] Configured coulomb_field_limit:"
        << "\n  raw internal value = "
        << configured_coulomb_field_limit
        << "\n  converted back to V/cm = "
        << Units::convert(
            configured_coulomb_field_limit,
            "V/cm"
        )
        << " V/cm"
        << "\n  stored coulomb_field_limit_ = "
        << coulomb_field_limit_
        << "\n  stored value converted back to V/cm = "
        << Units::convert(
            coulomb_field_limit_,
            "V/cm"
        )
        << " V/cm";
    // coulomb_field_limit_squared_ = config.get<double>("coulomb_field_limit") * config.get<double>("coulomb_field_limit") * 1e-10; // Convert from (V/cm)^2 to (MV/mm)^2 (internal field units)

    output_plots_ = config_.get<bool>("output_plots");
    output_linegraphs_ = config_.get<bool>("output_linegraphs");
    output_linegraphs_collected_ = config_.get<bool>("output_linegraphs_collected");
    output_linegraphs_recombined_ = config_.get<bool>("output_linegraphs_recombined");
    output_linegraphs_trapped_ = config_.get<bool>("output_linegraphs_trapped");
    output_rms_ = config_.get<bool>("output_rms");
    output_plots_step_ = config_.get<double>("output_plots_step");
    // Rickard 2026-04-05: Get step size for outputting propagation summary objects
    output_propagation_summary_step_ = config_.get<double>("output_propagation_summary_step");

    // Enable multithreading of this module if multithreading is enabled and no per-event output plots are requested:
    // FIXME: Review if this is really the case or we can still use multithreading
    if(!(config_.get<bool>("output_animations") || output_linegraphs_ || output_rms_)) {
        allow_multithreading();
    } else {
        LOG(WARNING) << "Per-event line graphs or animations requested, disabling parallel event processing";
    }

    // Parameter for charge transport in magnetic field (approximated from graphs:
    // http://www.ioffe.ru/SVA/NSM/Semicond/Si/electric.html) FIXME
    electron_Hall_ = 1.15;
    hole_Hall_ = 0.9;
}

// Copied from TransientPropagation.cpp with one section added to calculate the z-positions of the electrodes
void InteractivePropagationModule::initialize() {

    // Check for electric field
    if(!detector_->hasElectricField()) {
        LOG(WARNING) << "This detector does not have an electric field.";
    }

    if(!detector_->hasWeightingPotential()) {
        throw ModuleError("This module requires a weighting potential.");
    }

    if(detector_->getElectricFieldType() == FieldType::LINEAR) {
        LOG(ERROR) << "This module will likely produce unphysical results when applying linear electric fields.";
    }

    // Apply warnings if things are disabled
    if(!enable_diffusion_){
        LOG(WARNING) << "Diffusion is disabled in propagation. Results will be unphysical.";
    }

    if(!enable_coulomb_repulsion_){
        LOG(WARNING) << "Coulomb Repulsion has been disabled. Use TransientPropagation instead for this use case.";
    }

    // Prepare mobility model
    mobility_ = Mobility(config_, model_->getSensorMaterial(), detector_->hasDopingProfile());

    // Prepare recombination model
    recombination_ = Recombination(config_, detector_->hasDopingProfile());

    // Prepare trapping model
    trapping_ = Trapping(config_);

    // Prepare detrapping model
    detrapping_ = Detrapping(config_);

    // Impact ionization model
    multiplication_ = ImpactIonization(config_);

    // Calculate the locations of the upper and lower electrodes.
    // Assumes the local position origin is in the center of the detector
    if (include_mirror_charges_){

        // Determine the distance from the model origin to electrodes in both z-directions
        auto model_size = model_->getSize();
        LOG(DEBUG) << "Model size: " << model_size.x() << ", " << model_size.y() << ", " << model_size.z();
        
        auto model_center = model_->getModelCenter();
        LOG(DEBUG) << "Model center: " << model_center.x() << ", " << model_center.y() << ", " << model_center.z();

        auto sensor_center = model_->getSensorCenter();
        LOG(DEBUG) << "Sensor center: " << sensor_center.x() << ", " << sensor_center.y() << ", " << sensor_center.z();

        auto matrix_center = model_->getMatrixCenter();
        LOG(DEBUG) << "matrix center: " << matrix_center.x() << ", " << matrix_center.y() << ", " << matrix_center.z();

        z_lim_neg_ = model_center.z() - model_size.z() / 2;
        z_lim_pos_ = model_center.z() + model_size.z() / 2;

        //TODO: Correct algorithm to not assume it's in the center
            // z_lim_neg_ = model_->intersect(ROOT::Math::XYZVector(0, 0, -1), ROOT::Math::XYZPoint(0,0,0));
    }

    // Check multiplication and step size larger than a picosecond:
    if(!multiplication_.is<NoImpactIonization>() && timestep_ > 0.001) {
        LOG(WARNING) << "Charge multiplication enabled with maximum timestep larger than 1ps" << std::endl
                     << "This might lead to unphysical gain values.";
    }

    // Check for magnetic field
    has_magnetic_field_ = detector_->hasMagneticField();
    if(has_magnetic_field_) {
        if(config_.get<bool>("ignore_magnetic_field")) {
            has_magnetic_field_ = false;
            LOG(WARNING) << "A magnetic field is switched on, but is set to be ignored for this module.";
        } else {
            LOG(DEBUG) << "This detector sees a magnetic field.";
        }
    }

    if(output_plots_) {

        auto pitch_x = static_cast<double>(Units::convert(model_->getPixelSize().x(), "um"));
        auto pitch_y = static_cast<double>(Units::convert(model_->getPixelSize().y(), "um"));

        potential_difference_ = CreateHistogram<TH1D>(
            "potential_difference",
            "Weighting potential difference between two steps;#left|#Delta#phi_{w}#right| [a.u.];events",
            500,
            0,
            1);
        induced_charge_histo_ = CreateHistogram<TH1D>("induced_charge_histo",
                                                      "Induced charge per time, all pixels;Drift time [ns];charge [e]",
                                                      static_cast<int>(integration_time_ / timestep_),
                                                      0,
                                                      static_cast<double>(Units::convert(integration_time_, "ns")));
        induced_charge_e_histo_ =
            CreateHistogram<TH1D>("induced_charge_e_histo",
                                  "Induced charge per time, electrons only, all pixels;Drift time [ns];charge [e]",
                                  static_cast<int>(integration_time_ / timestep_),
                                  0,
                                  static_cast<double>(Units::convert(integration_time_, "ns")));
        induced_charge_h_histo_ =
            CreateHistogram<TH1D>("induced_charge_h_histo",
                                  "Induced charge per time, holes only, all pixels;Drift time [ns];charge [e]",
                                  static_cast<int>(integration_time_ / timestep_),
                                  0,
                                  static_cast<double>(Units::convert(integration_time_, "ns")));
        if(!multiplication_.is<NoImpactIonization>()) {
            induced_charge_primary_histo_ =
                CreateHistogram<TH1D>("induced_charge_primary_histo",
                                      "Induced charge per time, primaries only, all pixels;Drift time [ns];charge [e]",
                                      static_cast<int>(integration_time_ / timestep_),
                                      0,
                                      static_cast<double>(Units::convert(integration_time_, "ns")));
            induced_charge_primary_e_histo_ = CreateHistogram<TH1D>(
                "induced_charge_primary_e_histo",
                "Induced charge per time, primary electrons only, all pixels;Drift time [ns];charge [e]",
                static_cast<int>(integration_time_ / timestep_),
                0,
                static_cast<double>(Units::convert(integration_time_, "ns")));
            induced_charge_primary_h_histo_ =
                CreateHistogram<TH1D>("induced_charge_primary_h_histo",
                                      "Induced charge per time, primary holes only, all pixels;Drift time [ns];charge [e]",
                                      static_cast<int>(integration_time_ / timestep_),
                                      0,
                                      static_cast<double>(Units::convert(integration_time_, "ns")));
            induced_charge_secondary_histo_ =
                CreateHistogram<TH1D>("induced_charge_secondary_histo",
                                      "Induced charge per time, secondaries only, all pixels;Drift time [ns];charge [e]",
                                      static_cast<int>(integration_time_ / timestep_),
                                      0,
                                      static_cast<double>(Units::convert(integration_time_, "ns")));
            induced_charge_secondary_e_histo_ = CreateHistogram<TH1D>(
                "induced_charge_secondary_e_histo",
                "Induced charge per time, secondary electrons only, all pixels;Drift time [ns];charge [e]",
                static_cast<int>(integration_time_ / timestep_),
                0,
                static_cast<double>(Units::convert(integration_time_, "ns")));
            induced_charge_secondary_h_histo_ =
                CreateHistogram<TH1D>("induced_charge_secondary_h_histo",
                                      "Induced charge per time, secondary holes only, all pixels;Drift time [ns];charge [e]",
                                      static_cast<int>(integration_time_ / timestep_),
                                      0,
                                      static_cast<double>(Units::convert(integration_time_, "ns")));
        }
        induced_charge_vs_depth_histo_ =
            CreateHistogram<TH2D>("induced_charge_vs_depth_histo",
                                  "Induced charge per time vs depth, all pixels;Drift time [ns];depth [mm];charge [e]",
                                  static_cast<int>(integration_time_ / timestep_),
                                  0,
                                  static_cast<double>(Units::convert(integration_time_, "ns")),
                                  100,
                                  -model_->getSensorSize().z() / 2.,
                                  model_->getSensorSize().z() / 2.);
        induced_charge_e_vs_depth_histo_ = CreateHistogram<TH2D>(
            "induced_charge_e_vs_depth_histo",
            "Induced charge per time vs depth, electrons only, all pixels;Drift time [ns];depth [mm];charge [e]",
            static_cast<int>(integration_time_ / timestep_),
            0,
            static_cast<double>(Units::convert(integration_time_, "ns")),
            100,
            -model_->getSensorSize().z() / 2.,
            model_->getSensorSize().z() / 2.);
        induced_charge_h_vs_depth_histo_ = CreateHistogram<TH2D>(
            "induced_charge_h_vs_depth_histo",
            "Induced charge per time vs depth, holes only, all pixels;Drift time [ns];depth [mm];charge [e]",
            static_cast<int>(integration_time_ / timestep_),
            0,
            static_cast<double>(Units::convert(integration_time_, "ns")),
            100,
            -model_->getSensorSize().z() / 2.,
            model_->getSensorSize().z() / 2.);
        induced_charge_map_ = CreateHistogram<TH2D>(
            "induced_charge_map",
            "Induced charge as a function of in-pixel carrier position;x%pitch [#mum];y%pitch [#mum];charge [e]",
            static_cast<int>(pitch_x),
            -pitch_x / 2,
            pitch_x / 2,
            static_cast<int>(pitch_y),
            -pitch_y / 2,
            pitch_y / 2);
        induced_charge_e_map_ = CreateHistogram<TH2D>("induced_charge_e_map",
                                                      "Induced charge as a function of in-pixel carrier position, electrons "
                                                      "only;x%pitch [#mum];y%pitch [#mum];charge [e]",
                                                      static_cast<int>(pitch_x),
                                                      -pitch_x / 2,
                                                      pitch_x / 2,
                                                      static_cast<int>(pitch_y),
                                                      -pitch_y / 2,
                                                      pitch_y / 2);
        induced_charge_h_map_ = CreateHistogram<TH2D>(
            "induced_charge_h_map",
            "Induced charge as a function of in-pixel carrier position, holes only;x%pitch [#mum];y%pitch [#mum];charge [e]",
            static_cast<int>(pitch_x),
            -pitch_x / 2,
            pitch_x / 2,
            static_cast<int>(pitch_y),
            -pitch_y / 2,
            pitch_y / 2);

        step_length_histo_ =
            CreateHistogram<TH1D>("step_length_histo",
                                  "Step length;length [#mum];integration steps",
                                  100,
                                  0,
                                  static_cast<double>(Units::convert(0.25 * model_->getSensorSize().z(), "um")));
        group_size_histo_ = CreateHistogram<TH1D>("group_size_histo",
                                                  "Group size;size [charges];Number of groups",
                                                  static_cast<int>(100 * charge_per_step_),
                                                  0,
                                                  static_cast<int>(100 * charge_per_step_));

        drift_time_histo_ = CreateHistogram<TH1D>("drift_time_histo",
                                                  "Drift time;Drift time [ns];charge carriers",
                                                  static_cast<int>(Units::convert(integration_time_, "ns") * 5),
                                                  0,
                                                  static_cast<double>(Units::convert(integration_time_, "ns")));

        recombine_histo_ =
            CreateHistogram<TH1D>("recombination_histo",
                                  "Fraction of recombined charge carriers;recombination [N / N_{total}] ;number of events",
                                  100,
                                  0,
                                  1);
        recombination_time_histo_ =
            CreateHistogram<TH1D>("recombination_time_histo",
                                  "Time until recombination of charge carriers;time [ns];charge carriers",
                                  static_cast<int>(Units::convert(integration_time_, "ns") * 5),
                                  0,
                                  static_cast<double>(Units::convert(integration_time_, "ns")));
        trapped_histo_ = CreateHistogram<TH1D>(
            "trapping_histo", "Fraction of trapped charge carriers;trapping [N / N_{total}] ;number of events", 100, 0, 1);
        trapping_time_histo_ = CreateHistogram<TH1D>("trapping_time_histo",
                                                     "Local time of trapping of charge carriers;time [ns];charge carriers",
                                                     static_cast<int>(Units::convert(integration_time_, "ns") * 5),
                                                     0,
                                                     static_cast<double>(Units::convert(integration_time_, "ns")));
        detrapping_time_histo_ =
            CreateHistogram<TH1D>("detrapping_time_histo",
                                  "Time from trapping until detrapping of charge carriers;time [ns];charge carriers",
                                  static_cast<int>(Units::convert(integration_time_, "ns") * 5),
                                  0,
                                  static_cast<double>(Units::convert(integration_time_, "ns")));

        if(!multiplication_.is<NoImpactIonization>()) {
            gain_primary_histo_ = CreateHistogram<TH1D>(
                "gain_primary_histo",
                "Gain per primarily induced charge carrier group after propagation;gain;number of groups transported",
                24,
                1,
                25);
            gain_all_histo_ =
                CreateHistogram<TH1D>("gain_all_histo",
                                      "Gain per charge carrier group after propagation;gain;number of groups transported",
                                      24,
                                      1,
                                      25);
            gain_e_histo_ =
                CreateHistogram<TH1D>("gain_e_histo",
                                      "Gain per primary electron group after propagation;gain;number of groups transported",
                                      24,
                                      1,
                                      25);
            gain_h_histo_ =
                CreateHistogram<TH1D>("gain_h_histo",
                                      "Gain per primary hole group after propagation;gain;number of groups transported",
                                      24,
                                      1,
                                      25);
            multiplication_level_histo_ = CreateHistogram<TH1D>(
                "multiplication_level_histo",
                "Multiplication level of propagated charge carriers;multiplication level;charge carriers",
                static_cast<int>(max_multiplication_level_),
                0,
                static_cast<int>(max_multiplication_level_));
            multiplication_depth_histo_ =
                CreateHistogram<TH1D>("multiplication_depth_histo",
                                      "Generation depth of charge carriers via impact ionization;depth [mm];charge carriers",
                                      200,
                                      -model_->getSensorSize().z() / 2.,
                                      model_->getSensorSize().z() / 2.);
            gain_e_vs_x_ =
                CreateHistogram<TProfile>("gain_e_vs_x",
                                          "Gain per electron group after propagation vs x; x [mm]; gain per group",
                                          100,
                                          -model_->getSensorSize().x() / 2.,
                                          model_->getSensorSize().x() / 2.);
            gain_e_vs_y_ =
                CreateHistogram<TProfile>("gain_e_vs_y",
                                          "Gain per electron group after propagation vs y; x [mm]; gain per group",
                                          100,
                                          -model_->getSensorSize().y() / 2.,
                                          model_->getSensorSize().y() / 2.);
            gain_e_vs_z_ =
                CreateHistogram<TProfile>("gain_e_vs_z",
                                          "Gain per electron group after propagation vs z; x [mm]; gain per group",
                                          100,
                                          -model_->getSensorSize().z() / 2.,
                                          model_->getSensorSize().z() / 2.);
            gain_h_vs_x_ = CreateHistogram<TProfile>("gain_h_vs_x",
                                                     "Gain per hole group after propagation vs x; x [mm]; gain per group",
                                                     100,
                                                     -model_->getSensorSize().x() / 2.,
                                                     model_->getSensorSize().x() / 2.);
            gain_h_vs_y_ = CreateHistogram<TProfile>("gain_h_vs_y",
                                                     "Gain per hole group after propagation vs y; x [mm]; gain per group",
                                                     100,
                                                     -model_->getSensorSize().y() / 2.,
                                                     model_->getSensorSize().y() / 2.);
            gain_h_vs_z_ = CreateHistogram<TProfile>("gain_h_vs_z",
                                                     "Gain per hole group after propagation vs z; x [mm]; gain per group",
                                                     100,
                                                     -model_->getSensorSize().z() / 2.,
                                                     model_->getSensorSize().z() / 2.);
        }

        if (output_rms_){
            rms_total_graph_ = new TMultiGraph("rms_total_graph","Comparison of spread of electrons (dashed) and holes (solid);Drift time [ns];RMS [mm]");
            rms_e_subgraph_ = new TGraph();
            rms_e_subgraph_->SetNameTitle("rms_e_subgraph","Spread of electrons");
            rms_e_subgraph_->SetLineColor(kBlack);
            rms_e_subgraph_->SetLineStyle(kDashed);
            rms_h_subgraph_ = new TGraph();
            rms_h_subgraph_->SetNameTitle("rms_h_subgraph","Spread of holes");
            rms_h_subgraph_->SetLineColor(kBlack);
            rms_h_subgraph_->SetLineStyle(kSolid);

            rms_e_graph_ = new TMultiGraph("rms_e_graph","Spread of electrons(xyz=rgb);Drift time [ns];RMS [mm]");
            rms_x_e_subgraph_ = new TGraph();
            rms_x_e_subgraph_->SetNameTitle("rms_x_e_subgraph","Spread in X");
            rms_x_e_subgraph_->SetLineColor(kRed);
            rms_y_e_subgraph_ = new TGraph();
            rms_y_e_subgraph_->SetNameTitle("rms_y_e_subgraph","Spread in Y");
            rms_y_e_subgraph_->SetLineColor(kGreen);
            rms_z_e_subgraph_ = new TGraph();
            rms_z_e_subgraph_->SetNameTitle("rms_z_e_subgraph","Spread in Z");
            rms_z_e_subgraph_->SetLineColor(kBlue);
        }

        coulomb_mag_histo_ =
            CreateHistogram<TH1D>("coulomb_mag_histo",
                      "Direct Coulomb Field Interaction Magnitude;Interaction Field Magnitude [V/cm];Count",
                      200, // Number of bins for the field magnitude
                      0,   // Minimum field magnitude
                      Units::convert(
                        coulomb_field_limit_,
                        "V/cm"
                    )
            );
    }
}

// Copied from TransientPropagation.cpp with major modifications
void InteractivePropagationModule::run(Event* event) {
    auto deposits_message = messenger_->fetchMessage<DepositedChargeMessage>(this, event);

    // Create vector of propagated charges to output
    std::vector<PropagatedCharge> propagated_charges;

    // Create vector of propagation summaries to output
    std::vector<PropagationSummary> propagation_summaries;

    // List of points to plot to plot for output plots
    LineGraph::OutputPlotPoints output_plot_points;

    // Create vector of propagating charges to store each charge groups position, location, time, type, etc. at the start of propagation
    std::vector<PropagatedCharge> propagating_charges;

    //Create vector to temporarily store the applicable deposited charges
    std::vector<DepositedCharge> deposited_charges;

    //Storage of total charge
    unsigned int total_deposited_charge = 0;

    // Initial loop just to get the total charge of applicable charges
    for(const auto& deposit : deposits_message->getData()) {

        // Only process if within requested integration time: 
        if(deposit.getLocalTime() > integration_time_) {
            continue;
        }

        // Skip charges with type not included in propagation
        if(!propagate_electrons_ && deposit.getType() == allpix::CarrierType::ELECTRON){
            continue;
        }
        if(!propagate_holes_ && deposit.getType() == allpix::CarrierType::HOLE){
            continue;
        }

        total_deposits_++;
        total_deposited_charge += deposit.getCharge();

    }

    // The number of charges per charge group based on the total amount of charge
        // In the ideal case of a single deposit, this charge_per_step would give the desired number of charge groups
    auto default_charge_per_step = static_cast<unsigned int>(ceil(static_cast<double>(total_deposited_charge) / std::max(max_charge_groups_,static_cast<unsigned int>(1))));

    // Add the special 0 case where the max_charge_groups parameter is ignored
    if (max_charge_groups_ == 0) {
        default_charge_per_step = charge_per_step_;
    }

    // Choose the maximum charge_per_step value (which gives the minimum number of charge groups)
    auto charge_per_step = charge_per_step_;
    if (default_charge_per_step > charge_per_step) {
        charge_per_step = default_charge_per_step;
        LOG(INFO) << "max_charge_groups = " << max_charge_groups_ << " is the limiting factor in the charge_per_step = " << charge_per_step << " calculation";
    } else {
        LOG(INFO) << "charge_per_step = " << charge_per_step_ << " is the limiting factor in the charge_per_step = " << charge_per_step << " calculation (>" << default_charge_per_step << ")";
    }

    // Loop over all deposits for propagation
    for(const auto& deposit : deposits_message->getData()) {

        // Only process if within requested integration time: 
        if(deposit.getLocalTime() > integration_time_) {
            LOG(DEBUG) << "Skipping charge carriers deposited beyond integration time: "
                       << Units::display(deposit.getGlobalTime(), "ns") << " global / "
                       << Units::display(deposit.getLocalTime(), {"ns", "ps"}) << " local"
                       << "> Integration Time of " << integration_time_;
            continue;
        }

        // Skip charges with type not included in propagation
        if(!propagate_electrons_ && deposit.getType() == allpix::CarrierType::ELECTRON){
            LOG(DEBUG) << "Skipping "<< deposit.getCharge() <<" electron deposit as per configuration: " << Units::display(deposit.getLocalPosition(), {"mm", "um"});
            continue;
        }
        if(!propagate_holes_ && deposit.getType() == allpix::CarrierType::HOLE){
            LOG(DEBUG) << "Skipping "<< deposit.getCharge() <<" hole deposit as per configuration: " << Units::display(deposit.getLocalPosition(), {"mm", "um"});
            continue;
        }

        LOG(DEBUG) << "Set of "<<deposit.getCharge()<<" charge carriers (" << deposit.getType() << ") on "
                   << Units::display(deposit.getLocalPosition(), {"mm", "um"});
    
        // Split the deposit into charge groups
        unsigned int charges_remaining = deposit.getCharge();
        unsigned int charge_step = charge_per_step; // New charge step variable that gets reset for each deposit
        while(charges_remaining > 0) {

            // Define number of charges to be propagated and remove charges of this step from the total
            if(charge_per_step > charges_remaining) {
                charge_step = charges_remaining; // The final deposit has the remaining charge
            }
            charges_remaining -= charge_step;

            // Add charge to propagating charge vector to be time-stepped later
            PropagatedCharge propagating_charge(deposit.getLocalPosition(),
                                            deposit.getGlobalPosition(),
                                            deposit.getType(),
                                            charge_step,
                                            deposit.getLocalTime(), // The local deposition time
                                            deposit.getGlobalTime(), // The global deposition time
                                            CarrierState::MOTION,
                                            &deposit);

            propagating_charges.push_back(std::move(propagating_charge));
        }
    }
    
    if(max_charge_groups_ != 0 &&
    propagating_charges.size() > max_charge_groups_) {

        LOG(WARNING)
            << "Number of charge groups ("
            << propagating_charges.size()
            << ") exceeded set limit of "
            << max_charge_groups_
            << " due to the large number of deposits with low "
            "charge quantity (true limit = set limit + number "
            "of deposits)";
    }

    if(propagating_charges.empty()) {
        LOG(INFO)
            << "No applicable charge groups were produced for this event";
    } else {
        const double average_charges_per_group =
            static_cast<double>(
                total_deposited_charge
            )
            / static_cast<double>(
                propagating_charges.size()
            );

        LOG(INFO)
            << "Average number of charges per group is "
            << average_charges_per_group
            << " ("
            << propagating_charges.size()
            << " total)";
    }

    // Diagnose the temporal structure that a fully coupled RK4
    // implementation must handle. Splitting one deposit into many
    // charge groups produces repeated entries with the same activation
    // time, so the number of unique times can be much smaller than the
    // number of groups.
    std::map<double, unsigned int>
        groups_by_deposition_time;

    for(const auto& propagating_charge :
        propagating_charges) {

        groups_by_deposition_time[
            propagating_charge.getLocalTime()
        ]++;
    }

    if(!groups_by_deposition_time.empty()) {

        const double earliest_deposition_time =
            groups_by_deposition_time.begin()->first;

        const double latest_deposition_time =
            groups_by_deposition_time.rbegin()->first;

        unsigned int groups_active_at_zero = 0;
        unsigned int groups_active_by_midpoint = 0;
        unsigned int groups_active_before_step_end = 0;

        for(const auto& [
                deposition_time,
                number_of_groups
            ] : groups_by_deposition_time) {

            if(deposition_time <= 0.0) {
                groups_active_at_zero +=
                    number_of_groups;
            }

            if(deposition_time <=
            0.5 * timestep_) {

                groups_active_by_midpoint +=
                    number_of_groups;
            }

            if(deposition_time <
            timestep_) {

                groups_active_before_step_end +=
                    number_of_groups;
            }
        }

        const double deposition_time_span =
            latest_deposition_time
            - earliest_deposition_time;

        LOG(WARNING)
            << "[COUPLED_RK4_ACTIVATION_SUMMARY]"
            << "\n  total charge groups = "
            << propagating_charges.size()
            << "\n  unique exact deposition times = "
            << groups_by_deposition_time.size()
            << "\n  earliest deposition time = "
            << Units::convert(
                earliest_deposition_time,
                "ns"
            )
            << " ns"
            << "\n  latest deposition time = "
            << Units::convert(
                latest_deposition_time,
                "ns"
            )
            << " ns"
            << "\n  deposition-time span = "
            << Units::convert(
                deposition_time_span,
                "ns"
            )
            << " ns"
            << "\n  configured timestep = "
            << Units::convert(
                timestep_,
                "ns"
            )
            << " ns"
            << "\n  span / timestep = "
            << deposition_time_span
                / timestep_
            << "\n  groups active at t = 0 = "
            << groups_active_at_zero
            << "\n  groups active by K2/K3 midpoint = "
            << groups_active_by_midpoint
            << "\n  groups deposited before first step end = "
            << groups_active_before_step_end;

        unsigned int displayed_times = 0;

        for(const auto& [
                deposition_time,
                number_of_groups
            ] : groups_by_deposition_time) {

            if(displayed_times >= 20) {
                break;
            }

            LOG(WARNING)
                << "[COUPLED_RK4_ACTIVATION_TIME]"
                << "\n  activation index = "
                << displayed_times
                << "\n  deposition time = "
                << std::setprecision(
                    std::numeric_limits<double>::max_digits10
                )
                << Units::convert(
                    deposition_time,
                    "ns"
                )
                << " ns"
                << "\n  groups at this time = "
                << number_of_groups;

            displayed_times++;
        }
    }

    LOG(WARNING)
        << "[COULOMB_DEBUG_EVENT]"
        << "\n  event number = "
        << event->number
        << "\n  total deposited charge = "
        << total_deposited_charge
        << " e"
        << "\n  charge_per_step configured = "
        << charge_per_step_
        << "\n  charge_per_step actually used = "
        << charge_per_step
        << "\n  max_charge_groups = "
        << max_charge_groups_
        << "\n  number of propagating groups = "
        << propagating_charges.size()
        << "\n  Coulomb enabled = "
        << enable_coulomb_repulsion_
        << "\n  distance cutoff, internal = "
        << std::sqrt(
            coulomb_distance_limit_squared_
        )
        << "\n  distance cutoff = "
        << Units::convert(
            std::sqrt(
                coulomb_distance_limit_squared_
            ),
            "um"
        )
        << " um";
    // Propagation occurs within the following function call
    auto [recombined_charges_count, trapped_charges_count, propagated_charges_count] = propagate_together(event, propagating_charges, propagated_charges, propagation_summaries, output_plot_points);

    // Output plots if required
    if(output_linegraphs_) {
        LineGraph::Create(event->number, this, config_, output_plot_points, CarrierState::UNKNOWN);
        if(output_linegraphs_collected_) {
            LineGraph::Create(event->number, this, config_, output_plot_points, CarrierState::HALTED);
        }
        if(output_linegraphs_recombined_) {
            LineGraph::Create(event->number, this, config_, output_plot_points, CarrierState::RECOMBINED);
        }
        if(output_linegraphs_trapped_) {
            LineGraph::Create(event->number, this, config_, output_plot_points, CarrierState::TRAPPED);
        }
        if(config_.get<bool>("output_animations")) {
            LineGraph::Animate(event->number, this, config_, output_plot_points);
        }
    }

    LOG(INFO) << "Propagated " << propagated_charges_count << " charges" << std::endl
              << "Recombined " << recombined_charges_count << " charges during transport" << std::endl
              << "Trapped " << trapped_charges_count << " charges during transport";

    if(output_plots_) {
        auto total = (propagated_charges_count + recombined_charges_count + trapped_charges_count);
        recombine_histo_->Fill(static_cast<double>(recombined_charges_count) / (total == 0 ? 1 : total));
        trapped_histo_->Fill(static_cast<double>(trapped_charges_count) / (total == 0 ? 1 : total));
    }

    // Create a new message with propagated charges
    auto propagated_charge_message = std::make_shared<PropagatedChargeMessage>(std::move(propagated_charges), detector_);

    // Dispatch the message with propagated charges
    messenger_->dispatchMessage(this, std::move(propagated_charge_message), event);

    // Rickard 2026-04-05: Add message dispatching of propagation summaries for use in other modules or for output, if enabled in the configuration file will contain the summary for each sampling step
    if(output_propagation_summary_) {
        // Create a new message with propagation summaries
        auto propagation_summary_message =
            std::make_shared<PropagationSummaryMessage>(std::move(propagation_summaries), detector_);
        // Dispatch the message with propagation summaries
        messenger_->dispatchMessage(this, std::move(propagation_summary_message), event);
    }
}

// This function takes a list of propagating charges to propagate synchronously and places them in the propagated vector
std::tuple<unsigned int, unsigned int, unsigned int>
InteractivePropagationModule::propagate_together(Event* event,
                                                 std::vector<PropagatedCharge>& propagating_charges,
                                                 std::vector<PropagatedCharge>& propagated_charges,
                                                 std::vector<PropagationSummary>& propagation_summaries,
                                                 LineGraph::OutputPlotPoints& output_plot_points) const {

    unsigned int propagated_charges_count = 0;
    unsigned int recombined_charges_count = 0;
    unsigned int trapped_charges_count = 0;

    std::chrono::duration<double, std::nano> time_spent_coulomb{};

    // Define a function to compute the diffusion
    auto carrier_diffusion = [&](double efield_mag, double doping, double timestep, allpix::CarrierType type) -> Eigen::Vector3d {
        double diffusion_constant = boltzmann_kT_ * mobility_(type, efield_mag, doping);
        double diffusion_std_dev = std::sqrt(2. * diffusion_constant * timestep);

        // Compute the independent diffusion in three
        allpix::normal_distribution<double> gauss_distribution(0, diffusion_std_dev);
        auto x = gauss_distribution(event->getRandomEngine());
        auto y = gauss_distribution(event->getRandomEngine());
        auto z = gauss_distribution(event->getRandomEngine());
        return {x, y, z};
    };

    // Survival probability of this charge carrier package, evaluated at every step
    allpix::uniform_real_distribution<double> uniform_distribution(0, 1);

    // Create vectors that store charge information in a place that can be modified each time step 
        // They need to be here since they are used in dynamic field function, but they are set to initial states below
        // The order of objects within them must stay consistent
    std::vector<ROOT::Math::XYZPoint> charge_locations; // Current position of each charge 
    std::vector<ROOT::Math::XYZPoint> previous_charge_locations; // Positions of each charge at the previous time step (only updated once at the end of each timestep)
    std::vector<double> charge_times; // Most recent time for all of the charges (by the end of propagation they should all be aligned)
    std::vector<allpix::CarrierState> charge_states; // The state of propagation of each charge group (whether it's propagated, trapped, or halted)
    std::vector<allpix::CarrierState> previous_charge_states; // States frozen at the beginning of the current outer timestep. These are used for all Coulomb evaluations during that timestep.
    std::vector<std::uint8_t> all_sources_enabled;// Explicit source-activation mask. The existing independent solver enables every group here and continues to use deposition-time gating. The later coupled solver will provide a fixed mask for each substep.
    double_t time = 0; // The current time threshold (we only propagate charges near this time)

    // Temporary diagnostic counter:
    unsigned int coulomb_debug_call_count = 0; // Checking the calls

    // Temporary runtime diagnostics:
    bool coulomb_debug_first_pair_logged = false;
    bool coulomb_debug_first_cap_logged = false;
    enum class SourceActivationMode : std::uint8_t {
        DEPOSITION_TIME,
        EXPLICIT_MASK
    };

    // Numerical separation used only when two distinct charge groups occupy
    // exactly the same position. This preserves the previous regularization
    // distance while making the field evaluation deterministic.
    const double overlap_distance =
        std::sqrt(1e-15); // Internal distance unit: mm

    // SplitMix64 is used only to construct stable pair-dependent directions.
    // It does not modify the event random-number engine.
    const auto splitmix64 =
        [](std::uint64_t value) -> std::uint64_t {
            value += 0x9e3779b97f4a7c15ULL;
            value =
                (value ^ (value >> 30U))
                * 0xbf58476d1ce4e5b9ULL;
            value =
                (value ^ (value >> 27U))
                * 0x94d049bb133111ebULL;

            return value ^ (value >> 31U);
        };

    // Return a deterministic target-minus-source separation vector for an
    // exactly overlapping pair.
    //
    // The same unordered pair receives the same axis, while exchanging target
    // and source reverses the vector. This preserves pairwise antisymmetry.
    const auto deterministic_overlap_separation =
        [&](unsigned int target_index,
            unsigned int source_index)
            -> ROOT::Math::XYZVector {

        const auto lower_index =
            std::min(
                target_index,
                source_index
            );

        const auto upper_index =
            std::max(
                target_index,
                source_index
            );

        const auto lower_key =
            static_cast<std::uint64_t>(
                lower_index
            );

        const auto upper_key =
            static_cast<std::uint64_t>(
                upper_index
            );

        const auto event_key =
            static_cast<std::uint64_t>(
                event->number
            );

        std::uint64_t pair_key =
            (lower_key << 32U)
            | upper_key;

        pair_key ^=
            splitmix64(event_key);

        const auto hash_1 =
            splitmix64(pair_key);

        const auto hash_2 =
            splitmix64(hash_1);

        // Convert the upper 53 bits to values in [0, 1).
        constexpr double inverse_2_pow_53 =
            1.0 / 9007199254740992.0;

        const double uniform_1 =
            static_cast<double>(
                hash_1 >> 11U
            )
            * inverse_2_pow_53;

        const double uniform_2 =
            static_cast<double>(
                hash_2 >> 11U
            )
            * inverse_2_pow_53;

        // Uniform direction on the unit sphere:
        // z is uniform in [-1, 1], while phi is uniform in [0, 2pi).
        const double direction_z =
            2.0 * uniform_1 - 1.0;

        const double transverse_magnitude =
            std::sqrt(
                std::max(
                    0.0,
                    1.0
                        - direction_z
                        * direction_z
                )
            );

        const double phi =
            2.0
            * ROOT::Math::Pi()
            * uniform_2;

        // target - source must reverse when target and source are exchanged.
        const double orientation =
            target_index < source_index
                ? -1.0
                : 1.0;

        return ROOT::Math::XYZVector(
            orientation
                * overlap_distance
                * transverse_magnitude
                * std::cos(phi),
            orientation
                * overlap_distance
                * transverse_magnitude
                * std::sin(phi),
            orientation
                * overlap_distance
                * direction_z
        );
    };

    // Computes the coulomb force component of the e-field given a desired local point
    auto coulomb_efield =
        [&](double evaluation_time,
            const ROOT::Math::XYZPoint& point,
            unsigned int target_index,
            const std::vector<ROOT::Math::XYZPoint>& source_positions,
            const std::vector<allpix::CarrierState>& source_states,
            const std::vector<std::uint8_t>& source_active,
            SourceActivationMode activation_mode,
            bool record_diagnostics)
            -> Eigen::Vector3d {

        const bool debug_this_call =
            record_diagnostics
            && coulomb_debug_call_count < 20;

        unsigned int debug_sources_total = 0;
        unsigned int debug_sources_future = 0;
        unsigned int debug_sources_halted = 0;
        unsigned int debug_sources_recombined = 0;
        unsigned int debug_sources_self = 0;
        unsigned int debug_sources_zero_distance = 0;
        unsigned int debug_sources_outside_cutoff = 0;
        unsigned int debug_sources_eligible = 0;
        unsigned int debug_sources_overlapping = 0;
        unsigned int debug_sources_inactive = 0;


        auto coulomb_start = std::chrono::system_clock::now();

        ROOT::Math::XYZVector field = ROOT::Math::XYZVector(0,0,0);

        // Predefine some variables to eliminate redefinitions during loop
        ROOT::Math::XYZPoint local_position;
        int q;
        int sign;
        ROOT::Math::XYZVector dist_vector;
        double dist_mag2;
        double dist_mag;
        double interaction_magnitude;

        // Skip function entirely if disabled by configuration file
        // if (!enable_coulomb_repulsion_){
        //     return Eigen::Vector3d(field.x(),field.y(),field.z());
        // }
        if(!enable_coulomb_repulsion_) {
            if(debug_this_call) {
                LOG(WARNING)
                    << "[COULOMB_DEBUG_CALL]"
                    << "\n  call number = "
                    << coulomb_debug_call_count
                    << "\n  Coulomb repulsion is disabled";
            }

            if(record_diagnostics) {
                coulomb_debug_call_count++;
            }

            return Eigen::Vector3d(
                field.x(),
                field.y(),
                field.z()
            );
        }

        if(
            source_positions.size()
                != propagating_charges.size()
            || source_states.size()
                != propagating_charges.size()
            || source_active.size()
                != propagating_charges.size()
            || target_index
                >= propagating_charges.size()
        ) {
            throw ModuleError(
                "InteractivePropagation internal vector size or target-index "
                "mismatch in coulomb_efield"
            );
        }

        for(unsigned int i = 0; i < source_positions.size(); i++) {

            // TODO: Add check with (oc)tree object that only looks at charges within a certain distance

            // Only get fields from charges that have deposition time less than the current time (skip the ones that haven't been deposited yet)
            // This means that trapped charges at future times are okay, but not charges that haven't been deposited yet
            // Charges that have reached the sensor (HALTED) are assumed to be swept away and don't contribute to the coulomb field either
            // if (propagating_charges[i].getLocalTime() > time || charge_states[i] == allpix::CarrierState::HALTED || charge_states[i] == allpix::CarrierState::RECOMBINED){
            //     continue;
            // }
            debug_sources_total++;

            // Keep the source population fixed throughout one coupled
            // substep, including its K4 endpoint.
            if(source_active[i] == 0U) {
                debug_sources_inactive++;
                continue;
            }

            if(
                activation_mode
                    == SourceActivationMode::DEPOSITION_TIME
                && propagating_charges[i].getLocalTime()
                    > evaluation_time
            ) {
                debug_sources_future++;
                continue;
            }

            if(source_states[i] == CarrierState::HALTED) {
                debug_sources_halted++;
                continue;
            }

            if(source_states[i] == CarrierState::RECOMBINED) {
                debug_sources_recombined++;
                continue;
            }

            // Detect an exact overlap between two distinct charge groups.
            // Do not modify the physical source position here: mirror-charge
            // calculations should continue to use the actual source position.
            local_position =
                source_positions[i];

            const bool source_overlaps_target =
                target_index != i
                && local_position == point;

            if(source_overlaps_target) {
                debug_sources_overlapping++;
            }

            // Get the correct signed charge
            q = static_cast<int>(propagating_charges[i].getCharge()); // Positive charge [e]
            sign = static_cast<int8_t>(propagating_charges[i].getType()); // Sign of the charge q

            // Calculate the coulomb field due to charges that aren't the current charge
                // The calculation needs to be in the if-statement rather than a terminiation/continue since we still want to include the mirror charges of the current charge
            if(target_index == i) {
                debug_sources_self++;
            } else {

                if(source_overlaps_target) {
                    dist_vector =
                        deterministic_overlap_separation(
                            target_index,
                            i
                        );
                } else {
                    dist_vector =
                        point - local_position;
                }

                dist_mag2 =
                    dist_vector.Mag2();

                if(dist_mag2 <= 0.0) {
                    debug_sources_zero_distance++;
                } else if(
                    dist_mag2
                    >= coulomb_distance_limit_squared_
                ) {
                    debug_sources_outside_cutoff++;
                } else {

                    debug_sources_eligible++;

                    dist_mag =
                        ROOT::Math::sqrt(dist_mag2);

                    const auto uncapped_interaction_magnitude =
                        coulomb_K_
                        / relative_permittivity_
                        * static_cast<double>(q)
                        / dist_mag2;

                    interaction_magnitude =
                        std::min(
                            coulomb_field_limit_,
                            uncapped_interaction_magnitude
                        );

                    if(
                        record_diagnostics
                        && !coulomb_debug_first_pair_logged
                    ) {
                        LOG(WARNING)
                            << "[COULOMB_DEBUG_PAIR]"
                            << "\n  call number = "
                            << coulomb_debug_call_count
                            << "\n  field evaluation time = "
                            << Units::convert(
                                evaluation_time,
                                "ns"
                            )
                            << " ns"
                            << "\n  outer timestep start = "
                            << Units::convert(
                                time,
                                "ns"
                            )
                            << " ns"
                            << "\n  source group index = "
                            << i
                            << "\n  target group index = "
                            << target_index
                            << "\n  source charge = "
                            << q
                            << " e"
                            << "\n  separation = "
                            << Units::convert(
                                dist_mag,
                                "um"
                            )
                            << " um"
                            << "\n  separation squared = "
                            << dist_mag2
                            << " mm^2"
                            << "\n  cutoff = "
                            << Units::convert(
                                std::sqrt(
                                    coulomb_distance_limit_squared_
                                ),
                                "um"
                            )
                            << " um"
                            << "\n  relative permittivity = "
                            << relative_permittivity_
                            << "\n  uncapped field, internal = "
                            << uncapped_interaction_magnitude
                            << "\n  uncapped field = "
                            << Units::convert(
                                uncapped_interaction_magnitude,
                                "V/cm"
                            )
                            << " V/cm"
                            << "\n  stored cap, internal = "
                            << coulomb_field_limit_
                            << "\n  stored cap = "
                            << Units::convert(
                                coulomb_field_limit_,
                                "V/cm"
                            )
                            << " V/cm"
                            << "\n  field used, internal = "
                            << interaction_magnitude
                            << "\n  field used = "
                            << Units::convert(
                                interaction_magnitude,
                                "V/cm"
                            )
                            << " V/cm"
                            << "\n  cap applied = "
                            << (
                                uncapped_interaction_magnitude
                                    > coulomb_field_limit_
                                ? "YES"
                                : "NO"
                            );

                        coulomb_debug_first_pair_logged = true;
                    }

                    if(
                        record_diagnostics
                        && !coulomb_debug_first_cap_logged
                        && uncapped_interaction_magnitude
                        > coulomb_field_limit_
                    ) {
                        LOG(WARNING)
                            << "[COULOMB_DEBUG_CAP]"
                            << "\n  call number = "
                            << coulomb_debug_call_count
                            << "\n  source group index = "
                            << i
                            << "\n  target group index = "
                            << target_index
                            << "\n  separation = "
                            << Units::convert(
                                dist_mag,
                                "um"
                            )
                            << " um"
                            << "\n  source charge = "
                            << q
                            << " e"
                            << "\n  uncapped field = "
                            << Units::convert(
                                uncapped_interaction_magnitude,
                                "V/cm"
                            )
                            << " V/cm"
                            << "\n  cap = "
                            << Units::convert(
                                coulomb_field_limit_,
                                "V/cm"
                            )
                            << " V/cm"
                            << "\n  field used = "
                            << Units::convert(
                                interaction_magnitude,
                                "V/cm"
                            )
                            << " V/cm";

                        coulomb_debug_first_cap_logged = true;
                    }

                    if(
                        record_diagnostics
                        && output_plots_
                        && std::isfinite(
                            interaction_magnitude
                        )
                        && interaction_magnitude >= 0.0
                    ) {
                        coulomb_mag_histo_->Fill(
                            Units::convert(
                                interaction_magnitude,
                                "V/cm"
                            )
                        );
                    }

                    field =
                        field
                        + dist_vector
                            / dist_mag
                            * sign
                            * interaction_magnitude;
                }
            }

            // Skip mirror charges when specified
            if (!include_mirror_charges_){
                continue;
            }

            // Perform same for the mirror charges based on electrode positions (z_lim_neg and z_lim_pos)
                // Note: this assumes a parallel plate sensor (symmetry about z) in order to simplify the poisson equation to the mirror charge solution 
                // (potential is constant on each plate)
            ROOT::Math::XYZPoint mirror_position_neg = ROOT::Math::XYZPoint(local_position.x(), local_position.y(), 2*z_lim_neg_ - local_position.z());
            ROOT::Math::XYZPoint mirror_position_pos = ROOT::Math::XYZPoint(local_position.x(), local_position.y(), 2*z_lim_pos_ - local_position.z());

            // Apply field for negative-side mirror charge
            dist_vector = point - mirror_position_neg;
            dist_mag2 = dist_vector.Mag2();

            if(dist_mag2 > 0.0 && dist_mag2 < coulomb_distance_limit_squared_) {
                dist_mag = ROOT::Math::sqrt(dist_mag2);
                field = field - dist_vector / dist_mag * sign *
                        std::min(coulomb_field_limit_, coulomb_K_ / relative_permittivity_ * q / dist_mag2);
            }

            // Apply field for positive-side mirror charge
            dist_vector = point - mirror_position_pos;
            dist_mag2 = dist_vector.Mag2();

            if(dist_mag2 > 0.0 && dist_mag2 < coulomb_distance_limit_squared_) {
                dist_mag = ROOT::Math::sqrt(dist_mag2);
                field = field - dist_vector / dist_mag * sign * std::min(coulomb_field_limit_, coulomb_K_ / relative_permittivity_ * q / dist_mag2);
            }
        }

        // TODO: Rather than using coulomb_field_limit for each interaction, we could use it on the final value. 
            // (commented out for now since determining a good value is tricky)
            // if (field.mag2() > coulomb_field_limit_squared_){
            //     double field_mag = ROOT::Math::sqrt(field.mag2());
            //     // LOG(INFO) << "Skipping field with magnitude " << field_mag << " > "<<coulomb_field_limit_;
            //     field = Eigen::Vector3d(field.x() / field_mag * coulomb_field_limit_, field.y() / field_mag * coulomb_field_limit_, field.z() / field_mag * coulomb_field_limit_);
            //     // LOG(INFO) << "   now " << field.mag2();
            //     numSamePos+=1;
            // }

        if(debug_this_call) {
            LOG(WARNING)
                << "[COULOMB_DEBUG_CALL]"
                << "\n  call number = "
                << coulomb_debug_call_count
                << "\n  field evaluation time = "
                << Units::convert(
                    evaluation_time,
                    "ns"
                )
                << " ns"
                << "\n  outer timestep start = "
                << Units::convert(
                    time,
                    "ns"
                )
                << " ns"
                << "\n  target deposition time = "
                << Units::convert(
                    propagating_charges[target_index]
                        .getLocalTime(),
                    "ns"
                )
                << " ns"
                << "\n  target group index = "
                << target_index
                << "\n  number of propagating groups = "
                << propagating_charges.size()
                << "\n  Coulomb enabled = "
                << enable_coulomb_repulsion_
                << "\n  cutoff, internal = "
                << std::sqrt(
                    coulomb_distance_limit_squared_
                )
                << "\n  cutoff = "
                << Units::convert(
                    std::sqrt(
                        coulomb_distance_limit_squared_
                    ),
                    "um"
                )
                << " um"
                << "\n  sources inspected = "
                << debug_sources_total
                << "\n  skipped: inactive by explicit mask = "
                << debug_sources_inactive
                << "\n  skipped: future deposition = "
                << debug_sources_future
                << "\n  skipped: halted = "
                << debug_sources_halted
                << "\n  skipped: recombined = "
                << debug_sources_recombined
                << "\n  skipped: target self = "
                << debug_sources_self
                << "\n  exact zero distance = "
                << debug_sources_zero_distance
                << "\n  overlapping and offset = "
                << debug_sources_overlapping
                << "\n  outside cutoff = "
                << debug_sources_outside_cutoff
                << "\n  eligible direct interactions = "
                << debug_sources_eligible
                << "\n  resulting field vector, internal = ("
                << field.x()
                << ", "
                << field.y()
                << ", "
                << field.z()
                << ")"
                << "\n  resulting field vector in V/cm = ("
                << Units::convert(
                    field.x(),
                    "V/cm"
                )
                << ", "
                << Units::convert(
                    field.y(),
                    "V/cm"
                )
                << ", "
                << Units::convert(
                    field.z(),
                    "V/cm"
                )
                << ")"
                << "\n  resulting total field, internal = "
                << std::sqrt(
                    field.Mag2()
                )
                << "\n  resulting total field = "
                << Units::convert(
                    std::sqrt(
                        field.Mag2()
                    ),
                    "V/cm"
                )
                << " V/cm";
        }

        if(record_diagnostics) {
            coulomb_debug_call_count++;
        }

        Eigen::Vector3d output = Eigen::Vector3d(field.x(),field.y(),field.z());

        if(record_diagnostics) {
            auto coulomb_end = std::chrono::system_clock::now();
            time_spent_coulomb += coulomb_end - coulomb_start;
        }

        return output;
    };

    // Calculate the velocity of one target charge using an explicitly supplied
    // common source-position state.
    std::function<
        Eigen::Vector3d(
            double,
            const Eigen::Vector3d&,
            allpix::CarrierType,
            unsigned int,
            const std::vector<ROOT::Math::XYZPoint>&,
            const std::vector<allpix::CarrierState>&,
            const std::vector<std::uint8_t>&,
            SourceActivationMode,
            bool
        )
    > carrier_velocity_noB =
        [&](double evaluation_time,
            const Eigen::Vector3d& cur_pos,
            allpix::CarrierType type,
            unsigned int target_index,
            const std::vector<ROOT::Math::XYZPoint>& source_positions,
            const std::vector<allpix::CarrierState>& source_states,
            const std::vector<std::uint8_t>& source_active,
            SourceActivationMode activation_mode,
            bool record_coulomb_diagnostics)
            -> Eigen::Vector3d {

        const auto local_position =
            static_cast<ROOT::Math::XYZPoint>(
                cur_pos
            );

        const auto raw_field =
            detector_->getElectricField(
                local_position
            );

        Eigen::Vector3d efield(
            raw_field.x(),
            raw_field.y(),
            raw_field.z()
        );

        efield += coulomb_efield(
            evaluation_time,
            local_position,
            target_index,
            source_positions,
            source_states,
            source_active,
            activation_mode,
            record_coulomb_diagnostics
        );

        const auto doping =
            detector_->getDopingConcentration(
                local_position
            );

        return static_cast<int>(type)
            * mobility_(
                type,
                efield.norm(),
                doping
            )
            * efield;
    };

    std::function<
        Eigen::Vector3d(
            double,
            const Eigen::Vector3d&,
            allpix::CarrierType,
            unsigned int,
            const std::vector<ROOT::Math::XYZPoint>&,
            const std::vector<allpix::CarrierState>&,
            const std::vector<std::uint8_t>&,
            SourceActivationMode,
            bool
        )
    > carrier_velocity_withB =
        [&](double evaluation_time,
            const Eigen::Vector3d& cur_pos,
            allpix::CarrierType type,
            unsigned int target_index,
            const std::vector<ROOT::Math::XYZPoint>& source_positions,
            const std::vector<allpix::CarrierState>& source_states,
            const std::vector<std::uint8_t>& source_active,
            SourceActivationMode activation_mode,
            bool record_coulomb_diagnostics)
            -> Eigen::Vector3d {

        const auto local_position =
            static_cast<ROOT::Math::XYZPoint>(
                cur_pos
            );

        const auto raw_field =
            detector_->getElectricField(
                local_position
            );

        Eigen::Vector3d efield(
            raw_field.x(),
            raw_field.y(),
            raw_field.z()
        );

        efield += coulomb_efield(
            evaluation_time,
            local_position,
            target_index,
            source_positions,
            source_states,
            source_active,
            activation_mode,
            record_coulomb_diagnostics
        );

        const auto magnetic_field =
            detector_->getMagneticField(
                local_position
            );

        const Eigen::Vector3d bfield(
            magnetic_field.x(),
            magnetic_field.y(),
            magnetic_field.z()
        );

        const auto doping =
            detector_->getDopingConcentration(
                local_position
            );

        const auto mob =
            mobility_(
                type,
                efield.norm(),
                doping
            );

        const auto exb =
            efield.cross(
                bfield
            );

        const double hall_factor =
            type == CarrierType::ELECTRON
                ? electron_Hall_
                : hole_Hall_;

        const Eigen::Vector3d term1 =
            static_cast<int>(type)
            * mob
            * hall_factor
            * exb;

        const Eigen::Vector3d term2 =
            mob
            * mob
            * hall_factor
            * hall_factor
            * efield.dot(bfield)
            * bfield;

        const auto normalization =
            1.0
            + mob
            * mob
            * hall_factor
            * hall_factor
            * bfield.dot(bfield);

        return static_cast<int>(type)
            * mob
            * (
                efield
                + term1
                + term2
            )
            / normalization;
    };

    // Helper functions that convert between ROOT::Math:XYZPoint and Eigen::Vector3d (Eigen::Matrix<double, 3, 1>)
    auto convertPointToVector = [] (ROOT::Math::XYZPoint point) -> Eigen::Vector3d {
        return Eigen::Vector3d(point.x(), point.y(), point.z());
    };

    auto convertVectorToPoint = [] (Eigen::Vector3d vector) -> ROOT::Math::XYZPoint {
        return ROOT::Math::XYZPoint(vector.x(), vector.y(), vector.z());
    };

    auto convertRootVectorToEigenVector = [] (ROOT::Math::XYZVector vector) -> Eigen::Vector3d {
        return Eigen::Vector3d(vector.x(), vector.y(), vector.z());
    };

    // Create the pixel map which is used to collect the pulse objects
    std::vector< std::map<Pixel::Index, Pulse> > pixel_map_vector;

    // Create list of RK4 objects that correspond to each particle
    std::vector<allpix::RungeKutta<double, 4, 3>> runge_kutta_vector;

    // Initialize all vectors and the temporary independent RK4 objects.
    //
    // These RK4 objects still use the frozen beginning-of-step source state.
    // This preserves the existing behaviour while preparing the field and
    // velocity functions for the later globally coupled RK4 implementation.
    for(unsigned int i = 0;
        i < propagating_charges.size();
        i++) {

        const auto& charge =
            propagating_charges[i];

        const auto charge_type =
            charge.getType();

        std::function<
            Eigen::Matrix<double, 3, 1>(
                double,
                Eigen::Matrix<double, 3, 1>
            )
        > step_function;

        if(has_magnetic_field_) {
            step_function =
                [&, i, charge_type](
                    double evaluation_time,
                    Eigen::Vector3d trial_position
                ) -> Eigen::Vector3d {

                return carrier_velocity_withB(
                    evaluation_time,
                    trial_position,
                    charge_type,
                    i,
                    previous_charge_locations,
                    previous_charge_states,
                    all_sources_enabled,
                    SourceActivationMode::DEPOSITION_TIME,
                    true
                );
            };
        } else {
            step_function =
                [&, i, charge_type](
                    double evaluation_time,
                    Eigen::Vector3d trial_position
                ) -> Eigen::Vector3d {

                return carrier_velocity_noB(
                    evaluation_time,
                    trial_position,
                    charge_type,
                    i,
                    previous_charge_locations,
                    previous_charge_states,
                    all_sources_enabled,
                    SourceActivationMode::DEPOSITION_TIME,
                    true
                );
            };
        }

        auto rk =
            make_runge_kutta(
                tableau::RK4,
                step_function,
                timestep_,
                convertPointToVector(
                    charge.getLocalPosition()
                )
            );

        // Start propagation at the local deposition time.
        rk.advanceTime(
            charge.getLocalTime()
        );

        runge_kutta_vector.push_back(
            std::move(rk)
        );

        pixel_map_vector.emplace_back();

        charge_locations.push_back(
            charge.getLocalPosition()
        );

        previous_charge_locations.push_back(
            charge.getLocalPosition()
        );

        charge_times.push_back(
            charge.getLocalTime()
        );

        charge_states.push_back(
            charge.getState()
        );

        previous_charge_states.push_back(
            charge.getState()
        );

        all_sources_enabled.push_back(
            1U
        );

        if(output_linegraphs_) {
            output_plot_points.emplace_back(
                std::make_tuple(
                    charge.getGlobalTime(),
                    charge.getCharge(),
                    charge.getType(),
                    CarrierState::MOTION
                ),
                std::vector<ROOT::Math::XYZPoint>()
            );
        }
    }

    if(propagating_charges.empty()) {
        return std::make_tuple(0U, 0U, 0U);
    }

    // Set up variables that are changed each loop
    Eigen::Vector3d efield{};
    allpix::PropagatedCharge charge = propagating_charges[0];
    ROOT::Math::XYZPoint position{}; // = ROOT::Math::XYZPoint();
    ROOT::Math::XYZPoint previous_position{}; // = ROOT::Math::XYZPoint();
    allpix::CarrierType type = charge.getType();
    allpix::CarrierState state{};

    // Continue time propagation until the integration time has been reached
    for(time = 0; time < integration_time_; time += timestep_) { // time is the threshold value for each iteration

        // Based on the desired output_plots_step, display integration progress and calculate rms if desired
        if(std::fmod(time, output_plots_step_) < timestep_){
            // TODO: Change output_plots_step implementation to not depend on floating point errors.

            // This code could be useful in making the output_plots_step more consistent and not dependent on floating point errors
                // auto time_idx = static_cast<size_t>(runge_kutta_vector[i].getTime() / output_plots_step_);
                // while(next_idx <= time_idx) {
                    
                //     // output_plot_points.at(output_plot_index).second.push_back(static_cast<ROOT::Math::XYZPoint>(charge.getLocalPosition()));
                //     next_idx = output_plot_points.at(i).second.size();
                // }

            LOG(DEBUG) << "Time has reached " << time << "ns of " << integration_time_ << "ns";

            for (unsigned int i = 0; i < propagating_charges.size(); i++){
                // Add the current position to the linegraph associated with the current charge
                
            }

            // Get RMS of the charge distribution
            if (output_rms_){

                // Start by calculating the mean
                double x_mean_e = 0; double y_mean_e = 0; double z_mean_e = 0;
                double x_mean_h = 0; double y_mean_h = 0; double z_mean_h = 0;
                
                double num_e = 0; 
                double num_h = 0;
                for (unsigned int i = 0; i < charge_locations.size(); i++){
                    auto location = charge_locations[i];

                    // TODO: Think about whether there are certain states or time conditions we want to remove from RMS calc (ie RECOMBINED)

                    if (propagating_charges[i].getType() == allpix::CarrierType::ELECTRON){
                        num_e++;
                        x_mean_e += location.x();
                        y_mean_e += location.y();
                        z_mean_e += location.z();
                    }else{
                        num_h++;
                        x_mean_h += location.x();
                        y_mean_h += location.y();
                        z_mean_h += location.z();
                    }
                }

                if (num_e > 0){
                    x_mean_e = x_mean_e/num_e;
                    y_mean_e = y_mean_e/num_e;
                    z_mean_e = z_mean_e/num_e;
                }
                if (num_h > 0){
                    x_mean_h = x_mean_h/num_h;
                    y_mean_h = y_mean_h/num_h;
                    z_mean_h = z_mean_h/num_h;
                }
                
                // Now sum the square of the residuals (split up into x, y and z)
                double res2_x_e = 0; double res2_y_e = 0; double res2_z_e = 0;
                double res2_x_h = 0; double res2_y_h = 0; double res2_z_h = 0;

                for (unsigned int i = 0; i < charge_locations.size(); i++){

                    auto location = charge_locations[i];

                    if (propagating_charges[i].getType() == allpix::CarrierType::ELECTRON){
                        res2_x_e += (location.x() - x_mean_e)*(location.x() - x_mean_e);
                        res2_y_e += (location.y() - y_mean_e)*(location.y() - y_mean_e);
                        res2_z_e += (location.z() - z_mean_e)*(location.z() - z_mean_e);
                    }else{
                        res2_x_h += (location.x() - x_mean_h)*(location.x() - x_mean_h);
                        res2_y_h += (location.y() - y_mean_h)*(location.y() - y_mean_h);
                        res2_z_h += (location.z() - z_mean_h)*(location.z() - z_mean_h);
                    }
                }

                // Divide by the total number of charges of each type
                double rms_total_e = 0; double rms_x_e = 0; double rms_y_e = 0; double rms_z_e = 0;
                if (num_e > 0){
                    rms_total_e = sqrt((res2_x_e + res2_y_e + res2_z_e)/num_e);
                    rms_x_e = sqrt(res2_x_e/num_e);
                    rms_y_e = sqrt(res2_y_e/num_e);
                    rms_z_e = sqrt(res2_z_e/num_e);
                }
                double rms_total_h = 0; // double rms_x_h = 0; double rms_y_h = 0; double rms_z_h = 0;
                if (num_h > 0){
                    rms_total_h = sqrt((res2_x_h + res2_y_h + res2_z_h)/num_h);
                    // double rms_x_h = sqrt(res2_x_h/num_h); // Holes are less important, so ignore the separation of axes
                    // double rms_y_h = sqrt(res2_y_h/num_h);
                    // double rms_z_h = sqrt(res2_z_h/num_h);
                }

                // Add to ROOT Graphs
                rms_x_e_subgraph_->AddPoint(time, rms_x_e);
                rms_y_e_subgraph_->AddPoint(time, rms_y_e);
                rms_z_e_subgraph_->AddPoint(time, rms_z_e);
                rms_e_subgraph_->AddPoint(time, rms_total_e);
                rms_h_subgraph_->AddPoint(time, rms_total_h);

            }
        }

        if(output_propagation_summary_ && std::fmod(time, output_propagation_summary_step_) < timestep_) {
            
            const double nan = std::numeric_limits<double>::quiet_NaN();

            double sum_q_e = 0.0;
            double sum_x_e = 0.0;
            double sum_y_e = 0.0;
            double sum_z_e = 0.0;

            double sum_q_h = 0.0;
            double sum_x_h = 0.0;
            double sum_y_h = 0.0;
            double sum_z_h = 0.0;

            bool have_electrons = false;
            bool have_holes = false;

            double min_x_e = nan, max_x_e = nan;
            double min_y_e = nan, max_y_e = nan;
            double min_z_e = nan, max_z_e = nan;

            double min_x_h = nan, max_x_h = nan;
            double min_y_h = nan, max_y_h = nan;
            double min_z_h = nan, max_z_h = nan;

            for(unsigned int i = 0; i < charge_locations.size(); i++) {

                if(propagating_charges[i].getLocalTime() > time) {
                    continue;
                }

                if(charge_states[i] == CarrierState::RECOMBINED ||
                   charge_states[i] == CarrierState::TRAPPED) {
                    continue;
                }



                const double q = static_cast<double>(propagating_charges[i].getCharge());
                const auto& location = charge_locations[i];

                if(propagating_charges[i].getType() == CarrierType::ELECTRON) {
                    sum_q_e += q;
                    sum_x_e += q * location.x();
                    sum_y_e += q * location.y();
                    sum_z_e += q * location.z();

                    if(!have_electrons) {
                        min_x_e = max_x_e = location.x();
                        min_y_e = max_y_e = location.y();
                        min_z_e = max_z_e = location.z();
                        have_electrons = true;
                    } else {
                        min_x_e = std::min(min_x_e, location.x());
                        max_x_e = std::max(max_x_e, location.x());
                        min_y_e = std::min(min_y_e, location.y());
                        max_y_e = std::max(max_y_e, location.y());
                        min_z_e = std::min(min_z_e, location.z());
                        max_z_e = std::max(max_z_e, location.z());
                    }
                } else if(propagating_charges[i].getType() == CarrierType::HOLE) {
                    sum_q_h += q;
                    sum_x_h += q * location.x();
                    sum_y_h += q * location.y();
                    sum_z_h += q * location.z();

                    if(!have_holes) {
                        min_x_h = max_x_h = location.x();
                        min_y_h = max_y_h = location.y();
                        min_z_h = max_z_h = location.z();
                        have_holes = true;
                    } else {
                        min_x_h = std::min(min_x_h, location.x());
                        max_x_h = std::max(max_x_h, location.x());
                        min_y_h = std::min(min_y_h, location.y());
                        max_y_h = std::max(max_y_h, location.y());
                        min_z_h = std::min(min_z_h, location.z());
                        max_z_h = std::max(max_z_h, location.z());
                    }
                }
            }

            double mean_x_e = nan, mean_y_e = nan, mean_z_e = nan;
            double mean_x_h = nan, mean_y_h = nan, mean_z_h = nan;

            if(have_electrons && sum_q_e > 0.0) {
                mean_x_e = sum_x_e / sum_q_e;
                mean_y_e = sum_y_e / sum_q_e;
                mean_z_e = sum_z_e / sum_q_e;
            }

            if(have_holes && sum_q_h > 0.0) {
                mean_x_h = sum_x_h / sum_q_h;
                mean_y_h = sum_y_h / sum_q_h;
                mean_z_h = sum_z_h / sum_q_h;
            }

            double var_x_e = 0.0, var_y_e = 0.0, var_z_e = 0.0;
            double var_x_h = 0.0, var_y_h = 0.0, var_z_h = 0.0;

            for(unsigned int i = 0; i < charge_locations.size(); i++) {

                if(propagating_charges[i].getLocalTime() > time) {
                    continue;
                }

                if(charge_states[i] == CarrierState::RECOMBINED ||
                   charge_states[i] == CarrierState::TRAPPED) {
                        continue;
                    }

                const double q = static_cast<double>(propagating_charges[i].getCharge());
                const auto& location = charge_locations[i];

                if(propagating_charges[i].getType() == CarrierType::ELECTRON && have_electrons && sum_q_e > 0.0) {
                    const double dx = location.x() - mean_x_e;
                    const double dy = location.y() - mean_y_e;
                    const double dz = location.z() - mean_z_e;

                    var_x_e += q * dx * dx;
                    var_y_e += q * dy * dy;
                    var_z_e += q * dz * dz;
                } else if(propagating_charges[i].getType() == CarrierType::HOLE && have_holes && sum_q_h > 0.0) {
                    const double dx = location.x() - mean_x_h;
                    const double dy = location.y() - mean_y_h;
                    const double dz = location.z() - mean_z_h;

                    var_x_h += q * dx * dx;
                    var_y_h += q * dy * dy;
                    var_z_h += q * dz * dz;
                }
            }

            const double rms_x_e = (have_electrons && sum_q_e > 0.0) ? std::sqrt(var_x_e / sum_q_e) : nan;
            const double rms_y_e = (have_electrons && sum_q_e > 0.0) ? std::sqrt(var_y_e / sum_q_e) : nan;
            const double rms_z_e = (have_electrons && sum_q_e > 0.0) ? std::sqrt(var_z_e / sum_q_e) : nan;
            const double rms_e_e = (have_electrons && sum_q_e > 0.0) ? std::sqrt((var_x_e + var_y_e + var_z_e) / sum_q_e) : nan;

            const double rms_x_h = (have_holes && sum_q_h > 0.0) ? std::sqrt(var_x_h / sum_q_h) : nan;
            const double rms_y_h = (have_holes && sum_q_h > 0.0) ? std::sqrt(var_y_h / sum_q_h) : nan;
            const double rms_z_h = (have_holes && sum_q_h > 0.0) ? std::sqrt(var_z_h / sum_q_h) : nan;
            const double rms_e_h = (have_holes && sum_q_h > 0.0) ? std::sqrt((var_x_h + var_y_h + var_z_h) / sum_q_h) : nan;

            propagation_summaries.emplace_back(time,
                                            have_electrons,
                                            have_holes,
                                            mean_x_e,
                                            mean_y_e,
                                            mean_z_e,
                                            rms_x_e,
                                            rms_y_e,
                                            rms_z_e,
                                            rms_e_e,
                                            min_x_e,
                                            max_x_e,
                                            min_y_e,
                                            max_y_e,
                                            min_z_e,
                                            max_z_e,
                                            mean_x_h,
                                            mean_y_h,
                                            mean_z_h,
                                            rms_x_h,
                                            rms_y_h,
                                            rms_z_h,
                                            rms_e_h,
                                            min_x_h,
                                            max_x_h,
                                            min_y_h,
                                            max_y_h,
                                            min_z_h,
                                            max_z_h);
        }

        // Freeze the source state used by all Coulomb evaluations during
        // this outer timestep. This prevents the result from depending on
        // the order in which target groups are propagated.
        for(unsigned int i = 0;
            i < charge_locations.size();
            i++) {

            previous_charge_locations[i] =
                charge_locations[i];

            previous_charge_states[i] =
                charge_states[i];
        }

        if(time == 0.0) {

            const double step_start =
                time;

            const double step_end =
                time + timestep_;

            std::vector<double> exact_boundaries;
            exact_boundaries.reserve(
                propagating_charges.size() + 2
            );

            exact_boundaries.push_back(
                step_start
            );

            for(const auto& charge_group :
                propagating_charges) {

                const double activation_time =
                    charge_group.getLocalTime();

                if(
                    activation_time > step_start
                    && activation_time < step_end
                ) {
                    exact_boundaries.push_back(
                        activation_time
                    );
                }
            }

            exact_boundaries.push_back(
                step_end
            );

            std::sort(
                exact_boundaries.begin(),
                exact_boundaries.end()
            );

            exact_boundaries.erase(
                std::unique(
                    exact_boundaries.begin(),
                    exact_boundaries.end()
                ),
                exact_boundaries.end()
            );

            const double boundary_merge_tolerance =
                64.0
                * std::numeric_limits<double>::epsilon()
                * std::max(
                    {
                        1.0,
                        std::abs(step_start),
                        std::abs(step_end)
                    }
                );

            std::vector<double> merged_boundaries;
            merged_boundaries.reserve(
                exact_boundaries.size()
            );

            for(const double boundary :
                exact_boundaries) {

                if(
                    merged_boundaries.empty()
                    || std::abs(
                        boundary
                        - merged_boundaries.back()
                    ) > boundary_merge_tolerance
                ) {
                    merged_boundaries.push_back(
                        boundary
                    );
                }
            }

            LOG(WARNING)
                << "[COUPLED_RK4_SUBSTEP_SCHEDULE]"
                << "\n  exact boundary count = "
                << exact_boundaries.size()
                << "\n  exact substep count = "
                << exact_boundaries.size() - 1
                << "\n  merged boundary count = "
                << merged_boundaries.size()
                << "\n  merged substep count = "
                << merged_boundaries.size() - 1
                << "\n  merge tolerance = "
                << Units::convert(
                    boundary_merge_tolerance,
                    "ns"
                )
                << " ns";

            for(std::size_t substep_index = 0;
                substep_index + 1 < merged_boundaries.size();
                substep_index++) {

                const double substep_start =
                    merged_boundaries[substep_index];

                const double substep_end =
                    merged_boundaries[substep_index + 1];

                const double substep_size =
                    substep_end - substep_start;

                unsigned int active_groups_at_start = 0;
                unsigned int exactly_active_groups_at_start = 0;
                unsigned int newly_active_groups = 0;

                std::vector<std::uint8_t> substep_source_active(
                    propagating_charges.size(),
                    0U
                );

                for(unsigned int source_index = 0;
                    source_index < propagating_charges.size();
                    source_index++) {

                    const auto& charge_group =
                        propagating_charges[source_index];

                    const double activation_time =
                        charge_group.getLocalTime();

                    // The merged substep mask is authoritative for the
                    // future coupled RK4 calculation.
                    if(
                        activation_time
                        <= substep_start
                        + boundary_merge_tolerance
                    ) {
                        substep_source_active[source_index] =
                            1U;

                        active_groups_at_start++;
                    }

                    // Count the groups that would pass the old exact
                    // deposition-time condition at this retained boundary.
                    if(
                        activation_time
                        <= substep_start
                    ) {
                        exactly_active_groups_at_start++;
                    }

                    // Count groups represented by this merged boundary.
                    if(
                        std::abs(
                            activation_time
                            - substep_start
                        ) <= boundary_merge_tolerance
                    ) {
                        newly_active_groups++;
                    }
                }

                LOG(WARNING)
                    << "[COUPLED_RK4_SUBSTEP]"
                    << "\n  substep index = "
                    << substep_index
                    << "\n  start = "
                    << std::setprecision(
                        std::numeric_limits<double>::max_digits10
                    )
                    << Units::convert(
                        substep_start,
                        "ns"
                    )
                    << " ns"
                    << "\n  end = "
                    << Units::convert(
                        substep_end,
                        "ns"
                    )
                    << " ns"
                    << "\n  size = "
                    << Units::convert(
                        substep_size,
                        "ns"
                    )
                    << " ns"
                    << "\n  newly active groups = "
                    << newly_active_groups
                    << "\n  active groups during substep = "
                    << active_groups_at_start
                    << "\n  groups passing exact time gate = "
                    << exactly_active_groups_at_start;

                if(substep_index == 1) {

                    constexpr unsigned int diagnostic_target_index =
                        0U;

                    const auto& diagnostic_target_position =
                        previous_charge_locations[
                            diagnostic_target_index
                        ];

                    const auto explicit_mask_field =
                        coulomb_efield(
                            substep_start,
                            diagnostic_target_position,
                            diagnostic_target_index,
                            previous_charge_locations,
                            previous_charge_states,
                            substep_source_active,
                            SourceActivationMode::EXPLICIT_MASK,
                            false
                        );

                    const auto deposition_time_field =
                        coulomb_efield(
                            substep_start,
                            diagnostic_target_position,
                            diagnostic_target_index,
                            previous_charge_locations,
                            previous_charge_states,
                            substep_source_active,
                            SourceActivationMode::DEPOSITION_TIME,
                            false
                        );

                    const auto activation_mode_difference =
                        explicit_mask_field
                        - deposition_time_field;

                    LOG(WARNING)
                        << "[COUPLED_RK4_ACTIVATION_MODE_CHECK]"
                        << "\n  substep index = "
                        << substep_index
                        << "\n  evaluation time = "
                        << std::setprecision(
                            std::numeric_limits<double>::max_digits10
                        )
                        << Units::convert(
                            substep_start,
                            "ns"
                        )
                        << " ns"
                        << "\n  target group index = "
                        << diagnostic_target_index
                        << "\n  enabled by merged mask = "
                        << active_groups_at_start
                        << "\n  passing exact time gate = "
                        << exactly_active_groups_at_start
                        << "\n  explicit-mask field, internal = ("
                        << explicit_mask_field.x()
                        << ", "
                        << explicit_mask_field.y()
                        << ", "
                        << explicit_mask_field.z()
                        << ")"
                        << "\n  deposition-time field, internal = ("
                        << deposition_time_field.x()
                        << ", "
                        << deposition_time_field.y()
                        << ", "
                        << deposition_time_field.z()
                        << ")"
                        << "\n  field difference, internal = ("
                        << activation_mode_difference.x()
                        << ", "
                        << activation_mode_difference.y()
                        << ", "
                        << activation_mode_difference.z()
                        << ")"
                        << "\n  field-difference magnitude = "
                        << Units::convert(
                            activation_mode_difference.norm(),
                            "V/cm"
                        )
                        << " V/cm";
                }

                // Construct one diagnostic common K2 source state for the
                // final activation substep. This remains a shadow calculation:
                // no physical positions, states, RK objects, pulses, or random
                // numbers are modified.
                if(
                    substep_index + 2
                    == merged_boundaries.size()
                ) {

                    constexpr unsigned int
                        diagnostic_target_index = 0U;

                    if(
                        substep_source_active[
                            diagnostic_target_index
                        ] == 0U
                        || previous_charge_states[
                            diagnostic_target_index
                        ] != CarrierState::MOTION
                    ) {
                        throw ModuleError(
                            "Diagnostic coupled RK4 target is not active "
                            "and moving in the final activation substep"
                        );
                    }

                    std::vector<Eigen::Vector3d>
                        shadow_k1(
                            propagating_charges.size(),
                            Eigen::Vector3d::Zero()
                        );

                    std::vector<ROOT::Math::XYZPoint>
                        shadow_stage2_positions =
                            previous_charge_locations;

                    unsigned int
                        shadow_active_moving_groups = 0;

                    double
                        shadow_displacement_squared_sum =
                            0.0;

                    double
                        shadow_max_displacement =
                            0.0;

                    for(unsigned int source_index = 0;
                        source_index
                            < propagating_charges.size();
                        source_index++) {

                        if(
                            substep_source_active[
                                source_index
                            ] == 0U
                        ) {
                            continue;
                        }

                        if(
                            previous_charge_states[
                                source_index
                            ] != CarrierState::MOTION
                        ) {
                            continue;
                        }

                        const auto initial_position =
                            convertPointToVector(
                                previous_charge_locations[
                                    source_index
                                ]
                            );

                        const auto source_type =
                            propagating_charges[
                                source_index
                            ].getType();

                        if(has_magnetic_field_) {
                            shadow_k1[source_index] =
                                carrier_velocity_withB(
                                    substep_start,
                                    initial_position,
                                    source_type,
                                    source_index,
                                    previous_charge_locations,
                                    previous_charge_states,
                                    substep_source_active,
                                    SourceActivationMode::
                                        EXPLICIT_MASK,
                                    false
                                );
                        } else {
                            shadow_k1[source_index] =
                                carrier_velocity_noB(
                                    substep_start,
                                    initial_position,
                                    source_type,
                                    source_index,
                                    previous_charge_locations,
                                    previous_charge_states,
                                    substep_source_active,
                                    SourceActivationMode::
                                        EXPLICIT_MASK,
                                    false
                                );
                        }

                        const auto stage2_position =
                            initial_position
                            + 0.5
                                * substep_size
                                * shadow_k1[
                                    source_index
                                ];

                        shadow_stage2_positions[
                            source_index
                        ] =
                            convertVectorToPoint(
                                stage2_position
                            );

                        const double displacement =
                            (
                                stage2_position
                                - initial_position
                            ).norm();

                        shadow_displacement_squared_sum +=
                            displacement
                            * displacement;

                        shadow_max_displacement =
                            std::max(
                                shadow_max_displacement,
                                displacement
                            );

                        shadow_active_moving_groups++;
                    }

                    const double
                        shadow_rms_displacement =
                            shadow_active_moving_groups > 0
                                ? std::sqrt(
                                    shadow_displacement_squared_sum
                                    / static_cast<double>(
                                        shadow_active_moving_groups
                                    )
                                )
                                : 0.0;

                    const auto&
                        diagnostic_target_stage2_position =
                            shadow_stage2_positions[
                                diagnostic_target_index
                            ];

                    // Both calls use the same target trial position and
                    // explicit activation mask. They differ only in whether
                    // the source cloud is frozen or moved to the common K2
                    // stage positions.
                    const auto
                        frozen_source_stage2_field =
                            coulomb_efield(
                                substep_start,
                                diagnostic_target_stage2_position,
                                diagnostic_target_index,
                                previous_charge_locations,
                                previous_charge_states,
                                substep_source_active,
                                SourceActivationMode::
                                    EXPLICIT_MASK,
                                false
                            );

                    const auto
                        common_source_stage2_field =
                            coulomb_efield(
                                substep_start,
                                diagnostic_target_stage2_position,
                                diagnostic_target_index,
                                shadow_stage2_positions,
                                previous_charge_states,
                                substep_source_active,
                                SourceActivationMode::
                                    EXPLICIT_MASK,
                                false
                            );

                    const auto
                        stage2_source_coupling_difference =
                            common_source_stage2_field
                            - frozen_source_stage2_field;

                    LOG(WARNING)
                        << "[COUPLED_RK4_COMMON_STAGE_CHECK]"
                        << "\n  substep index = "
                        << substep_index
                        << "\n  substep start = "
                        << std::setprecision(
                            std::numeric_limits<double>::
                                max_digits10
                        )
                        << Units::convert(
                            substep_start,
                            "ns"
                        )
                        << " ns"
                        << "\n  substep size = "
                        << Units::convert(
                            substep_size,
                            "ns"
                        )
                        << " ns"
                        << "\n  active moving groups = "
                        << shadow_active_moving_groups
                        << "\n  diagnostic target index = "
                        << diagnostic_target_index
                        << "\n  RMS K2 source displacement = "
                        << Units::convert(
                            shadow_rms_displacement,
                            "um"
                        )
                        << " um"
                        << "\n  maximum K2 source displacement = "
                        << Units::convert(
                            shadow_max_displacement,
                            "um"
                        )
                        << " um"
                        << "\n  frozen-source K2 field, internal = ("
                        << frozen_source_stage2_field.x()
                        << ", "
                        << frozen_source_stage2_field.y()
                        << ", "
                        << frozen_source_stage2_field.z()
                        << ")"
                        << "\n  common-source K2 field, internal = ("
                        << common_source_stage2_field.x()
                        << ", "
                        << common_source_stage2_field.y()
                        << ", "
                        << common_source_stage2_field.z()
                        << ")"
                        << "\n  common-minus-frozen field, internal = ("
                        << stage2_source_coupling_difference.x()
                        << ", "
                        << stage2_source_coupling_difference.y()
                        << ", "
                        << stage2_source_coupling_difference.z()
                        << ")"
                        << "\n  field-difference magnitude = "
                        << Units::convert(
                            stage2_source_coupling_difference.norm(),
                            "V/cm"
                        )
                        << " V/cm";
                }
            }
        }

        // Move all charges by a single timestep
        for (unsigned int i = 0; i < propagating_charges.size(); i++){
            
            // Update local variables for convenient access and reduced array calling
            auto &runge_kutta = runge_kutta_vector[i];
            position = convertVectorToPoint(runge_kutta.getValue());
            state = charge_states[i];

            // TODO: Change output_plots_step implementation to not depend on floating point errors.
            if(output_linegraphs_ && std::fmod(time, output_plots_step_) < timestep_) {
                output_plot_points.at(i).second.push_back(position);
            }

            // Only propagate within a timestep range above the time threshold (time <= rk_time < time + timestep_)
            if(runge_kutta.getTime() < time || runge_kutta.getTime() >= time + timestep_){ 
                continue;
            }
            // Now the propagations are calculated only for those in the proper range

            if(state == CarrierState::TRAPPED){ 
                // If it reaches here, it must be within the time range and previously set to trapped. So, we can remove the trapped state and continue propagation
                state = CarrierState::MOTION;
            }else if(state == CarrierState::RECOMBINED || state == CarrierState::HALTED || state == CarrierState::UNKNOWN){
                // I don't think any charges will reach here since they would have to be advanced forward with one of these states.
                continue;
            }
            // At this point, the state must be MOTION and we continue with the propagation

            // Update more local variables that aren't needed above (saves this for after the time and state filtering)
            charge = propagating_charges[i];
            previous_position = previous_charge_locations[i];
            type = charge.getType();


            // Get electric field at current (pre-step) position
            // TODO: add a storage of the dynamic field so that we don't have to calculate it an extra time for use in diffusion
            const double pre_step_time =
                runge_kutta.getTime();

            efield =
                convertRootVectorToEigenVector(
                    detector_->getElectricField(
                        position
                    )
                );

            efield += coulomb_efield(
                pre_step_time,
                position,
                i,
                previous_charge_locations,
                previous_charge_states,
                all_sources_enabled,
                SourceActivationMode::DEPOSITION_TIME,
                true
            );
            auto doping = detector_->getDopingConcentration(position); //TODO: Does doping affect the dynamic field at all?

            // Execute a Runge-Kutta step and verify how the solver advances
            // its internal physical time.
            const double rk_time_before =
                runge_kutta.getTime();

            auto step =
                runge_kutta.step();

            const double rk_time_after =
                runge_kutta.getTime();

            if(time == 0.0 && i < 4) {
                LOG(WARNING)
                    << "[RK4_TIME_ADVANCE_DEBUG]"
                    << "\n  target group index = "
                    << i
                    << "\n  configured timestep = "
                    << Units::convert(
                        timestep_,
                        "ns"
                    )
                    << " ns"
                    << "\n  RK time before step = "
                    << Units::convert(
                        rk_time_before,
                        "ns"
                    )
                    << " ns"
                    << "\n  RK time after step = "
                    << Units::convert(
                        rk_time_after,
                        "ns"
                    )
                    << " ns"
                    << "\n  RK time increment = "
                    << Units::convert(
                        rk_time_after
                            - rk_time_before,
                        "ns"
                    )
                    << " ns";
            }

            charge_times[i] =
                rk_time_after;
            
            // Get the new position due to the electric field
            position = convertVectorToPoint(runge_kutta.getValue());

            // Apply diffusion step (if enabled)
            if (enable_diffusion_){
                auto diffusion = carrier_diffusion(efield.norm(), doping, timestep_, charge.getType());
                position = ROOT::Math::XYZPoint(position.x() + diffusion.x(), position.y() + diffusion.y(), position.z() + diffusion.z());
            }

            // If charge carrier reaches implant, interpolate surface position for higher accuracy:
            if(auto implant = model_->isWithinImplant(position)) {
                LOG(TRACE) << "Carrier in implant: " << Units::display(position, {"nm"});
                auto new_position = model_->getImplantIntercept(implant.value(), previous_position, position);
                position = convertPointToVector(new_position);
                state = CarrierState::HALTED;
                // The runge kutta's time will remain at the time that this gets triggered
            }

            // Check for overshooting outside the sensor and correct for it:
            if(!model_->isWithinSensor(position)) {
                // Reflect off the sensor surface with a certain probability, otherwise halt motion:
                if(uniform_distribution(event->getRandomEngine()) > surface_reflectivity_) {
                    LOG(TRACE) << "Carrier outside sensor: "
                            << Units::display(position, {"nm"});
                    state = CarrierState::HALTED;
                }

                auto intercept = model_->getSensorIntercept(previous_position, position);

                if(state == CarrierState::HALTED) {
                    position = intercept;
                } else {
                    // geom. reflection on x-y plane at upper sensor boundary (we have an implant on the lower edge)
                    position = ROOT::Math::XYZPoint(position.x(), position.y(), 2. * intercept.z() - position.z());
                    LOG(TRACE) << "Carrier was reflected on the sensor surface to "
                            << Units::display(position, {"um", "nm"});

                    // Re-check if we ended in an implant - corner case.
                    if(model_->isWithinImplant(position)) {
                        LOG(TRACE) << "Ended in implant after reflection - halting";
                        state = CarrierState::HALTED;
                    }

                    // Re-check if we are within the sensor - reflection at sensor side walls:
                    if(!model_->isWithinSensor(position)) {
                        position = intercept;
                        state = CarrierState::HALTED;
                    }
                }
                LOG(TRACE) << "Moved carrier to: " << Units::display(position, {"nm"});
            }

            // Update final position after applying corrections from surface intercepts
            runge_kutta.setValue(convertPointToVector(position));

            // Update position vector after e-field and diffusion so it is up to date in in dynamic field calculation
            charge_locations[i] = position;

            // Update step length histogram
            if(output_plots_) {
                step_length_histo_->Fill(static_cast<double>(Units::convert(step.value.norm(), "um")));
            }

            // Physics effects:

            // Apply stochastic bulk processes only while the carrier is
            // still moving. A collected or recombined carrier must not
            // subsequently become trapped.
            if(
                state == CarrierState::MOTION
                && recombination_(
                    type,
                    doping,
                    uniform_distribution(
                        event->getRandomEngine()
                    ),
                    timestep_
                )
            ) {
                state =
                    CarrierState::RECOMBINED;
            }

            if(
                state == CarrierState::MOTION
                && trapping_(
                    type,
                    uniform_distribution(
                        event->getRandomEngine()
                    ),
                    timestep_,
                    efield.norm()
                )
            ) {
                state =
                    CarrierState::TRAPPED;

                if(output_plots_) {
                    trapping_time_histo_->Fill(
                        runge_kutta.getTime(),
                        charge.getCharge()
                    );
                }

                const auto detrap_time =
                    detrapping_(
                        type,
                        uniform_distribution(
                            event->getRandomEngine()
                        ),
                        efield.norm()
                    );

                runge_kutta.advanceTime(
                    detrap_time
                );

                if(runge_kutta.getTime() <
                integration_time_) {

                    LOG(TRACE)
                        << "Charge carrier will detrap after "
                        << Units::display(
                            detrap_time,
                            {"ns", "us"}
                        );

                    if(output_plots_) {
                        detrapping_time_histo_->Fill(
                            static_cast<double>(
                                Units::convert(
                                    detrap_time,
                                    "ns"
                                )
                            ),
                            charge.getCharge()
                        );
                    }
                }
            }

            // No multiplication occurs since adding more charge groups increases simulation time dramatically

            // Signal calculation:

            // Find the nearest pixel - before and after the step
            auto [xpixel, ypixel] = model_->getPixelIndex(position);
            auto [last_xpixel, last_ypixel] = model_->getPixelIndex(previous_position);
            auto idx = Pixel::Index(xpixel, ypixel);
            auto neighbors = model_->getNeighbors(idx, distance_);

            // If the charge carrier crossed pixel boundaries, ensure that we always calculate the induced current for both of
            // them by extending the induction matrix temporarily. Otherwise we end up doing "double-counting" because we would
            // only jump "into" a pixel but never "out". At the border of the induction matrix, this would create an imbalance.
            if(last_xpixel != xpixel || last_ypixel != ypixel) {
                auto last_idx = Pixel::Index(last_xpixel, last_ypixel);
                neighbors.merge(model_->getNeighbors(last_idx, distance_));
                LOG(TRACE) << "Carrier crossed boundary from pixel " << Pixel::Index(last_xpixel, last_ypixel) << " to pixel "
                        << Pixel::Index(xpixel, ypixel);
            }
            LOG(TRACE) << "Moving carriers below pixel " << Pixel::Index(xpixel, ypixel) << " from "
                    << Units::display(previous_position, {"um", "mm"}) << " to "
                    << Units::display(position, {"um", "mm"}) << ", "
                    << Units::display(runge_kutta.getTime(), "ns");

            for(const auto& pixel_index : neighbors) {
                auto ramo = detector_->getWeightingPotential(position, pixel_index);
                auto last_ramo = detector_->getWeightingPotential(previous_position, pixel_index);

                // Induced charge on electrode is q_int = q * (phi(x1) - phi(x0))
                auto induced = charge.getCharge() * (ramo - last_ramo) * static_cast<std::underlying_type<CarrierType>::type>(type);

                // This line is commented since we are not applying multiplication
                // auto induced_primary = level != 0 ? 0.
                //                                 : initial_charge * (ramo - last_ramo) *
                //                                         static_cast<std::underlying_type<CarrierType>::type>(type);
                auto induced_primary = induced;
                auto induced_secondary = induced - induced_primary; // TODO: If multiplocation isn't reimplemented, remove the redundant info

                LOG(TRACE) << "Pixel " << pixel_index << " dPhi = " << (ramo - last_ramo) << ", induced " << type
                        << " q = " << Units::display(induced, "e");

                // Create pulse if it doesn't exist. Store induced charge in the returned pulse iterator
                auto pixel_map_iterator = pixel_map_vector[i].emplace(pixel_index, Pulse(timestep_, integration_time_));
                try {
                    pixel_map_iterator.first->second.addCharge(induced, runge_kutta.getTime());
                } catch(const PulseBadAllocException& e) {
                    LOG(ERROR) << e.what() << std::endl
                            << "Ignoring pulse contribution at time "
                            << Units::display(runge_kutta.getTime(), {"ms", "us", "ns"});
                }

                if(output_plots_) {
                    auto inPixel_um_x = (position.x() - model_->getPixelCenter(xpixel, ypixel).x()) * 1e3;
                    auto inPixel_um_y = (position.y() - model_->getPixelCenter(xpixel, ypixel).y()) * 1e3;

                    potential_difference_->Fill(std::fabs(ramo - last_ramo));
                    induced_charge_histo_->Fill(runge_kutta.getTime(), induced);
                    induced_charge_vs_depth_histo_->Fill(runge_kutta.getTime(), position.z(), induced);
                    induced_charge_map_->Fill(inPixel_um_x, inPixel_um_y, induced);
                    if(type == CarrierType::ELECTRON) {
                        induced_charge_e_histo_->Fill(runge_kutta.getTime(), induced);
                        induced_charge_e_vs_depth_histo_->Fill(
                            runge_kutta.getTime(), position.z(), induced);
                        induced_charge_e_map_->Fill(inPixel_um_x, inPixel_um_y, induced);
                    } else {
                        induced_charge_h_histo_->Fill(runge_kutta.getTime(), induced);
                        induced_charge_h_vs_depth_histo_->Fill(
                            runge_kutta.getTime(), position.z(), induced);
                        induced_charge_h_map_->Fill(inPixel_um_x, inPixel_um_y, induced);
                    }
                    if(!multiplication_.is<NoImpactIonization>()) { //TODO: If muliplication isn't reimplemented, remove the primary and secondary histogram
                        induced_charge_primary_histo_->Fill(runge_kutta.getTime(), induced_primary);
                        induced_charge_secondary_histo_->Fill(runge_kutta.getTime(), induced_secondary);
                        if(type == CarrierType::ELECTRON) {
                            induced_charge_primary_e_histo_->Fill(runge_kutta.getTime(), induced_primary);
                            induced_charge_secondary_e_histo_->Fill(runge_kutta.getTime(),
                                                                    induced_secondary);
                        } else {
                            induced_charge_primary_h_histo_->Fill(runge_kutta.getTime(), induced_primary);
                            induced_charge_secondary_h_histo_->Fill(runge_kutta.getTime(),
                                                                    induced_secondary);
                        }
                    }
                }
            }
            // Increase charge at the end of the step in case of impact ionization (commented since we are not performing multiplication)
                // charge += n_secondaries;

            // Set the values in vectors to keep them in sync with the propagation
            charge_states[i] = state;

        }

    }

    // Add final charges to propagated charges vector
    LOG(INFO) << "Outputing propagated charges";
    for (unsigned int i = 0; i < propagating_charges.size(); i++){

        charge = propagating_charges[i];
        auto runge_kutta = runge_kutta_vector[i];

        if(output_linegraphs_) {
            std::get<3>(output_plot_points.at(i).first) = charge_states[i];
        }

        // Create PropagatedCharge object and add it to the list
        auto local_position = convertVectorToPoint(runge_kutta.getValue());
        auto global_position = detector_->getGlobalPosition(local_position);
        auto local_time = runge_kutta.getTime();
        auto global_time = local_time - charge.getLocalTime() + charge.getGlobalTime();

        const DepositedCharge* deposit = charge.getDepositedCharge();

        PropagatedCharge propagated_charge(local_position,
                                        global_position,
                                        charge.getType(),
                                        std::move(pixel_map_vector[i]),
                                        local_time,
                                        global_time,
                                        charge_states[i],
                                        deposit);

        LOG(DEBUG) << " Propagated " << charge << " (initial: " << charge.getCharge() << ") to "
                << Units::display(local_position, {"mm", "um"}) << " in " << Units::display(runge_kutta.getTime(), "ns")
                << " time, induced " << Units::display(propagated_charge.getCharge(), {"e"})
                << ", final state: " << allpix::to_string(charge_states[i]);

        propagated_charges.push_back(std::move(propagated_charge));

        // Calculate the final totals for the recombined, trapped, and propagated charges
        
        if(charge_states[i] == CarrierState::RECOMBINED) {
            recombined_charges_count += charge.getCharge();
            if(output_plots_) {
                recombination_time_histo_->Fill(runge_kutta.getTime(), charge.getCharge());
            }
        } else if(charge_states[i] == CarrierState::TRAPPED) { // If the charge still has the TRAPPED state at the integration time, it is clear that the detrapping time was sufficiently large
            trapped_charges_count += charge.getCharge();
        } else {
            propagated_charges_count += charge.getCharge();
        }
    
        if(output_plots_) {
            drift_time_histo_->Fill(static_cast<double>(Units::convert(runge_kutta.getTime(), "ns")), charge.getCharge()); //TODO: Check whether we need to remove the "dead time" before deposition
            group_size_histo_->Fill(charge.getCharge());
        }
    }

    LOG(INFO) << "The running of the coulomb_efield function took a combined " << time_spent_coulomb.count()/1e6 << "ms";

    return std::make_tuple(recombined_charges_count,trapped_charges_count,propagated_charges_count);
}

// Copied from TransientPropagation.cpp with addition of rms plots
void InteractivePropagationModule::finalize() {
    LOG(INFO) << deposits_exceeding_max_groups_ * 100.0 / total_deposits_ << "% of deposits have charge exceeding the "
              << max_charge_groups_ << " charge groups allowed, with a charge_per_step value of " << charge_per_step_ << "."; //TODO: Change to make sense with the new interpretation of max_charge_groups
    if(output_plots_) {
        group_size_histo_->Get()->GetXaxis()->SetRange(1, group_size_histo_->Get()->GetNbinsX() + 1);

        potential_difference_->Write();
        step_length_histo_->Write();
        group_size_histo_->Write();
        drift_time_histo_->Write();
        recombine_histo_->Write();
        recombination_time_histo_->Write();
        trapped_histo_->Write();
        induced_charge_histo_->Write();
        induced_charge_e_histo_->Write();
        induced_charge_h_histo_->Write();
        if(!multiplication_.is<NoImpactIonization>()) {
            induced_charge_primary_histo_->Write();
            induced_charge_primary_e_histo_->Write();
            induced_charge_primary_h_histo_->Write();
            induced_charge_secondary_histo_->Write();
            induced_charge_secondary_e_histo_->Write();
            induced_charge_secondary_h_histo_->Write();
        }
        induced_charge_vs_depth_histo_->Write();
        induced_charge_e_vs_depth_histo_->Write();
        induced_charge_h_vs_depth_histo_->Write();
        induced_charge_map_->Write();
        induced_charge_e_map_->Write();
        induced_charge_h_map_->Write();
        if(!multiplication_.is<NoImpactIonization>()) {
            gain_primary_histo_->Write();
            gain_all_histo_->Write();
            gain_e_histo_->Write();
            gain_h_histo_->Write();
            multiplication_level_histo_->Write();
            multiplication_depth_histo_->Write();
            gain_e_vs_x_->Write();
            gain_e_vs_y_->Write();
            gain_e_vs_z_->Write();
            gain_h_vs_x_->Write();
            gain_h_vs_y_->Write();
            gain_h_vs_z_->Write();
        }

        if (output_rms_){
            rms_total_graph_->Add(rms_e_subgraph_);
            rms_total_graph_->Add(rms_h_subgraph_);

            rms_e_graph_->Add(rms_x_e_subgraph_);
            rms_e_graph_->Add(rms_y_e_subgraph_);
            rms_e_graph_->Add(rms_z_e_subgraph_);
            rms_e_graph_->Add(rms_e_subgraph_);

            rms_total_graph_->Write();
            rms_e_graph_->Write();
        }

        coulomb_mag_histo_->Write();
    }
}
