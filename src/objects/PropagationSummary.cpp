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
                                       double mean_x,
                                       double mean_y,
                                       double mean_z,
                                       double rms_x,
                                       double rms_y,
                                       double rms_z,
                                       double rms_e,
                                       double min_x,
                                       double max_x,
                                       double min_y,
                                       double max_y,
                                       double min_z,
                                       double max_z)
    : local_time_(local_time), mean_x_(mean_x), mean_y_(mean_y), mean_z_(mean_z), rms_x_(rms_x), rms_y_(rms_y),
      rms_z_(rms_z), rms_e_(rms_e), min_x_(min_x), max_x_(max_x), min_y_(min_y), max_y_(max_y), min_z_(min_z), max_z_(max_z) {}

double PropagationSummary::getLocalTime() const {
    return local_time_;
}

double PropagationSummary::getMeanX() const {
    return mean_x_;
}
double PropagationSummary::getMeanY() const {
    return mean_y_;
}
double PropagationSummary::getMeanZ() const {
    return mean_z_;
}

double PropagationSummary::getRMSX() const {
    return rms_x_;
}
double PropagationSummary::getRMSY() const {
    return rms_y_;
}
double PropagationSummary::getRMSZ() const {
    return rms_z_;
}
double PropagationSummary::getRMSE() const {
    return rms_e_;
}

double PropagationSummary::getMinX() const {
    return min_x_;
}

double PropagationSummary::getMaxX() const {
    return max_x_;
}

double PropagationSummary::getMinY() const {
    return min_y_;
}

double PropagationSummary::getMaxY() const {
    return max_y_;
}

double PropagationSummary::getMinZ() const {
    return min_z_;
}

double PropagationSummary::getMaxZ() const {
    return max_z_;
}

void PropagationSummary::print(std::ostream& out) const {
    out << "PropagationSummary:"
        << " t=" << local_time_ << " mean=(" << mean_x_ << ", " << mean_y_ << ", " << mean_z_ << ")"
        << " rms=(" << rms_x_ << ", " << rms_y_ << ", " << rms_z_ << ", " << rms_e_ << ")"
        << " bounds=(" << min_x_ << ", " << max_x_ << "), (" << min_y_ << ", " << max_y_ << "), (" << min_z_ << ", " << max_z_ << ")";
}

void PropagationSummary::loadHistory() {}

void PropagationSummary::petrifyHistory() {}

ClassImp(PropagationSummary)