/**
 * @file
 * @brief Implements a user hook for Geant4 stepping action
 *
 * @copyright Copyright (c) 2023-2025 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "StepInfoUserHookG4.hpp"

#include <memory>
#include <string>
#include <utility>

#include <G4Gamma.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "DepositionGeant4Module.hpp"

#include "core/messenger/Messenger.hpp"
#include "core/module/Event.hpp"
#include "core/module/Module.hpp"
#include "core/utils/log.h"
#include "objects/RayleighScatter.hpp"
#include "tools/geant4/RunManager.hpp"

using namespace allpix;

thread_local std::vector<RayleighScatter> StepInfoUserHookG4::rayleigh_scatters_;

void StepInfoUserHookG4::UserSteppingAction(const G4Step* aStep) {
    if(aStep->GetStepLength() < 0) {
        LOG(WARNING) << "Negative step length found; aborting event.";

        // Load the G4 run manager, which is owned by the geometry builder, and abort the current event:
        G4RunManager::GetRunManager()->AbortRun();
        return;
    }

    auto* track = aStep->GetTrack();
    if(track == nullptr || track->GetDefinition() != G4Gamma::GammaDefinition()) {
        return;
    }

    auto* pre_step = aStep->GetPreStepPoint();
    auto* post_step = aStep->GetPostStepPoint();

    if(pre_step == nullptr || post_step == nullptr) {
        return;
    }

    const auto* process = post_step->GetProcessDefinedStep();
    if(process == nullptr) {
        return;
    }

    const auto process_name = process->GetProcessName();

    // Geant4 process name for photon Rayleigh scattering is normally "Rayl".
    if(process_name != "Rayl") {
        return;
    }

    const auto position = post_step->GetPosition();

    // Momentum vectors before and after the Rayleigh scattering process.
    // Stored in MeV, with XYZ components in global coordinates.
    const auto incoming_momentum = pre_step->GetMomentum();
    const auto outgoing_momentum = post_step->GetMomentum();

    rayleigh_scatters_.emplace_back(
        ROOT::Math::XYZPoint(position.x() / mm, position.y() / mm, position.z() / mm),
        ROOT::Math::XYZVector(incoming_momentum.x() / MeV, incoming_momentum.y() / MeV, incoming_momentum.z() / MeV),
        ROOT::Math::XYZVector(outgoing_momentum.x() / MeV, outgoing_momentum.y() / MeV, outgoing_momentum.z() / MeV),
        post_step->GetGlobalTime() / ns,
        track->GetTrackID(),
        track->GetParentID());

    LOG(TRACE) << "Stored Rayleigh scatter for G4 track " << track->GetTrackID() << " at (" << position.x() / mm << ", "
               << position.y() / mm << ", " << position.z() / mm << ") mm";
}

void StepInfoUserHookG4::clearRayleighScatterInfo() {
    rayleigh_scatters_.clear();
}

std::size_t StepInfoUserHookG4::getStoredRayleighScatters() {
    return rayleigh_scatters_.size();
}

void StepInfoUserHookG4::dispatchRayleighScatterMessage(Module* module, Messenger* messenger, Event* event) {
    if(rayleigh_scatters_.empty()) {
        return;
    }

    LOG(DEBUG) << "Dispatching " << rayleigh_scatters_.size() << " RayleighScatter object(s)";

    auto rayleigh_scatter_message = std::make_shared<RayleighScatterMessage>(std::move(rayleigh_scatters_));
    messenger->dispatchMessage(module, std::move(rayleigh_scatter_message), event);

    rayleigh_scatters_.clear();
}