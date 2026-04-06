/**
 * @file
 * @brief Implementation of propagation summary object
 *
 * @copyright Copyright (c) 2017-2025 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "PropagationSummary.hpp"

using namespace allpix;

PropagationSummary::PropagationSummary(double local_time,
                                       bool has_electrons,
                                       bool has_holes,
                                       double mean_x_e,
                                       double mean_y_e,
                                       double mean_z_e,
                                       double rms_x_e,
                                       double rms_y_e,
                                       double rms_z_e,
                                       double rms_e_e,
                                       double min_x_e,
                                       double max_x_e,
                                       double min_y_e,
                                       double max_y_e,
                                       double min_z_e,
                                       double max_z_e,
                                       double mean_x_h,
                                       double mean_y_h,
                                       double mean_z_h,
                                       double rms_x_h,
                                       double rms_y_h,
                                       double rms_z_h,
                                       double rms_e_h,
                                       double min_x_h,
                                       double max_x_h,
                                       double min_y_h,
                                       double max_y_h,
                                       double min_z_h,
                                       double max_z_h)
    : local_time_(local_time), has_electrons_(has_electrons), has_holes_(has_holes), mean_x_e_(mean_x_e),
      mean_y_e_(mean_y_e), mean_z_e_(mean_z_e), rms_x_e_(rms_x_e), rms_y_e_(rms_y_e), rms_z_e_(rms_z_e),
      rms_e_e_(rms_e_e), min_x_e_(min_x_e), max_x_e_(max_x_e), min_y_e_(min_y_e), max_y_e_(max_y_e),
      min_z_e_(min_z_e), max_z_e_(max_z_e), mean_x_h_(mean_x_h), mean_y_h_(mean_y_h), mean_z_h_(mean_z_h),
      rms_x_h_(rms_x_h), rms_y_h_(rms_y_h), rms_z_h_(rms_z_h), rms_e_h_(rms_e_h), min_x_h_(min_x_h),
      max_x_h_(max_x_h), min_y_h_(min_y_h), max_y_h_(max_y_h), min_z_h_(min_z_h), max_z_h_(max_z_h) {}

double PropagationSummary::getLocalTime() const {
    return local_time_;
}

bool PropagationSummary::hasElectrons() const {
    return has_electrons_;
}

bool PropagationSummary::hasHoles() const {
    return has_holes_;
}

double PropagationSummary::getMeanXE() const {
    return mean_x_e_;
}
double PropagationSummary::getMeanYE() const {
    return mean_y_e_;
}
double PropagationSummary::getMeanZE() const {
    return mean_z_e_;
}

double PropagationSummary::getRMSXE() const {
    return rms_x_e_;
}
double PropagationSummary::getRMSYE() const {
    return rms_y_e_;
}
double PropagationSummary::getRMSZE() const {
    return rms_z_e_;
}
double PropagationSummary::getRMSEE() const {
    return rms_e_e_;
}

double PropagationSummary::getMinXE() const {
    return min_x_e_;
}
double PropagationSummary::getMaxXE() const {
    return max_x_e_;
}
double PropagationSummary::getMinYE() const {
    return min_y_e_;
}
double PropagationSummary::getMaxYE() const {
    return max_y_e_;
}
double PropagationSummary::getMinZE() const {
    return min_z_e_;
}
double PropagationSummary::getMaxZE() const {
    return max_z_e_;
}

double PropagationSummary::getMeanXH() const {
    return mean_x_h_;
}
double PropagationSummary::getMeanYH() const {
    return mean_y_h_;
}
double PropagationSummary::getMeanZH() const {
    return mean_z_h_;
}

double PropagationSummary::getRMSXH() const {
    return rms_x_h_;
}
double PropagationSummary::getRMSYH() const {
    return rms_y_h_;
}
double PropagationSummary::getRMSZH() const {
    return rms_z_h_;
}
double PropagationSummary::getRMSEH() const {
    return rms_e_h_;
}

double PropagationSummary::getMinXH() const {
    return min_x_h_;
}
double PropagationSummary::getMaxXH() const {
    return max_x_h_;
}
double PropagationSummary::getMinYH() const {
    return min_y_h_;
}
double PropagationSummary::getMaxYH() const {
    return max_y_h_;
}
double PropagationSummary::getMinZH() const {
    return min_z_h_;
}
double PropagationSummary::getMaxZH() const {
    return max_z_h_;
}

void PropagationSummary::print(std::ostream& out) const {
    out << "PropagationSummary:"
        << " t=" << local_time_
        << " has_electrons=" << has_electrons_
        << " has_holes=" << has_holes_
        << " mean_e=(" << mean_x_e_ << ", " << mean_y_e_ << ", " << mean_z_e_ << ")"
        << " rms_e=(" << rms_x_e_ << ", " << rms_y_e_ << ", " << rms_z_e_ << ", " << rms_e_e_ << ")"
        << " bounds_e=(" << min_x_e_ << ", " << max_x_e_ << "), (" << min_y_e_ << ", " << max_y_e_ << "), ("
        << min_z_e_ << ", " << max_z_e_ << ")"
        << " mean_h=(" << mean_x_h_ << ", " << mean_y_h_ << ", " << mean_z_h_ << ")"
        << " rms_h=(" << rms_x_h_ << ", " << rms_y_h_ << ", " << rms_z_h_ << ", " << rms_e_h_ << ")"
        << " bounds_h=(" << min_x_h_ << ", " << max_x_h_ << "), (" << min_y_h_ << ", " << max_y_h_ << "), ("
        << min_z_h_ << ", " << max_z_h_ << ")";
}

void PropagationSummary::loadHistory() {}

void PropagationSummary::petrifyHistory() {}

ClassImp(PropagationSummary)