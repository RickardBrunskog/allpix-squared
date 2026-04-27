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
                            double max_z_h);

        double getLocalTime() const;

        bool hasElectrons() const;
        bool hasHoles() const;

        double getMeanXE() const;
        double getMeanYE() const;
        double getMeanZE() const;

        double getRMSXE() const;
        double getRMSYE() const;
        double getRMSZE() const;
        double getRMSEE() const;

        double getMinXE() const;
        double getMaxXE() const;
        double getMinYE() const;
        double getMaxYE() const;
        double getMinZE() const;
        double getMaxZE() const;

        double getMeanXH() const;
        double getMeanYH() const;
        double getMeanZH() const;

        double getRMSXH() const;
        double getRMSYH() const;
        double getRMSZH() const;
        double getRMSEH() const;

        double getMinXH() const;
        double getMaxXH() const;
        double getMinYH() const;
        double getMaxYH() const;
        double getMinZH() const;
        double getMaxZH() const;

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
        ClassDefOverride(PropagationSummary, 2); // NOLINT

    private:
        double local_time_{};

        bool has_electrons_{false};
        bool has_holes_{false};

        double mean_x_e_{};
        double mean_y_e_{};
        double mean_z_e_{};

        double rms_x_e_{};
        double rms_y_e_{};
        double rms_z_e_{};
        double rms_e_e_{};

        double min_x_e_{};
        double max_x_e_{};
        double min_y_e_{};
        double max_y_e_{};
        double min_z_e_{};
        double max_z_e_{};

        double mean_x_h_{};
        double mean_y_h_{};
        double mean_z_h_{};

        double rms_x_h_{};
        double rms_y_h_{};
        double rms_z_h_{};
        double rms_e_h_{};

        double min_x_h_{};
        double max_x_h_{};
        double min_y_h_{};
        double max_y_h_{};
        double min_z_h_{};
        double max_z_h_{};
    };

    /**
     * @brief Typedef for message carrying propagation summaries
     */
    using PropagationSummaryMessage = Message<PropagationSummary>;

} // namespace allpix

#endif /* ALLPIX_PROPAGATION_SUMMARY_H */