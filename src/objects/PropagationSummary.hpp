/**
 * @file
 * @brief Definition of propagation summary object
 *
 * @copyright Copyright (c) 2017-2025 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef ALLPIX_PROPAGATION_SUMMARY_H
#define ALLPIX_PROPAGATION_SUMMARY_H

#include "Object.hpp"

namespace allpix {

    /**
     * @ingroup Objects
     * @brief Summary of propagated charge cloud properties at one sampled time point
     */
    class PropagationSummary : public Object {
    public:
        /**
         * @brief Default constructor for ROOT I/O
         */
        PropagationSummary() = default;

        /**
         * @brief Construct a propagation summary point
         * @param local_time Local propagation time of this summary sample
         * @param mean_x Charge-weighted mean x position
         * @param mean_y Charge-weighted mean y position
         * @param mean_z Charge-weighted mean z position
         * @param rms_x Charge-weighted RMS spread in x
         * @param rms_y Charge-weighted RMS spread in y
         * @param rms_z Charge-weighted RMS spread in z
         * @param rms_e Total 3D RMS spread
         * @param min_x Minimum x position of the propagated charge cloud
         * @param max_x Maximum x position of the propagated charge cloud
         * @param min_y Minimum y position of the propagated charge cloud
         * @param max_y Maximum y position of the propagated charge cloud
         * @param min_z Minimum z position of the propagated charge cloud
         * @param max_z Maximum z position of the propagated charge cloud
         */
        PropagationSummary(double local_time,
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
                           double max_z);

        double getLocalTime() const;

        double getMeanX() const;
        double getMeanY() const;
        double getMeanZ() const;

        double getRMSX() const;
        double getRMSY() const;
        double getRMSZ() const;
        double getRMSE() const;

        double getMinX() const;
        double getMaxX() const;
        double getMinY() const;
        double getMaxY() const;
        double getMinZ() const;
        double getMaxZ() const;

        /**
         * @brief Print an ASCII representation of PropagationSummary to the given stream
         * @param out Stream to print to
         */
        void print(std::ostream& out) const override;

        void loadHistory() override;
        void petrifyHistory() override;

        /**
         * @brief ROOT class definition
         */
        ClassDefOverride(PropagationSummary, 1); // NOLINT

    private:
        double local_time_{};
        double mean_x_{};
        double mean_y_{};
        double mean_z_{};

        double rms_x_{};
        double rms_y_{};
        double rms_z_{};
        double rms_e_{};
        double min_x_{};
        double max_x_{};
        double min_y_{};
        double max_y_{};
        double min_z_{};
        double max_z_{};
    };

    /**
     * @brief Typedef for message carrying propagation summaries
     */
    using PropagationSummaryMessage = Message<PropagationSummary>;

} // namespace allpix

#endif /* ALLPIX_PROPAGATION_SUMMARY_H */