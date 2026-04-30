/**
 * @file
 * @brief Implementation of Rayleigh scatter object
 *
 * @copyright Copyright (c) 2017-2025 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "RayleighScatter.hpp"

#include <ostream>
#include <utility>

using namespace allpix;

RayleighScatter::RayleighScatter(ROOT::Math::XYZPoint position,
                                 ROOT::Math::XYZVector incoming_momentum,
                                 ROOT::Math::XYZVector outgoing_momentum,
                                 double global_time,
                                 int track_id,
                                 int parent_id)
    : position_(std::move(position)), incoming_momentum_(std::move(incoming_momentum)),
      outgoing_momentum_(std::move(outgoing_momentum)), global_time_(global_time), track_id_(track_id),
      parent_id_(parent_id) {}

ROOT::Math::XYZPoint RayleighScatter::getPosition() const {
    return position_;
}

ROOT::Math::XYZVector RayleighScatter::getIncomingMomentum() const {
    return incoming_momentum_;
}

ROOT::Math::XYZVector RayleighScatter::getOutgoingMomentum() const {
    return outgoing_momentum_;
}

double RayleighScatter::getGlobalTime() const {
    return global_time_;
}

int RayleighScatter::getTrackID() const {
    return track_id_;
}

int RayleighScatter::getParentID() const {
    return parent_id_;
}

void RayleighScatter::print(std::ostream& out) const {
    out << "RayleighScatter:"
        << " track_id=" << track_id_
        << " parent_id=" << parent_id_
        << " time=" << global_time_ << " ns"
        << " position=(" << position_.X() << ", " << position_.Y() << ", " << position_.Z() << ") mm"
        << " incoming_momentum=(" << incoming_momentum_.X() << ", " << incoming_momentum_.Y() << ", "
        << incoming_momentum_.Z() << ") MeV"
        << " outgoing_momentum=(" << outgoing_momentum_.X() << ", " << outgoing_momentum_.Y() << ", "
        << outgoing_momentum_.Z() << ") MeV";
}

void RayleighScatter::loadHistory() {}

void RayleighScatter::petrifyHistory() {}

ClassImp(RayleighScatter)
