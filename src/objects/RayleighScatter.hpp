/**
 * @file
 * @brief Definition of Rayleigh scatter object
 *
 * @copyright Copyright (c) 2017-2025 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef ALLPIX_RAYLEIGH_SCATTER_H
#define ALLPIX_RAYLEIGH_SCATTER_H

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "Object.hpp"

namespace allpix {

    /**
     * @ingroup Objects
     * @brief Monte Carlo truth object storing one Rayleigh scattering interaction of a photon
     */
    class RayleighScatter : public Object {
    public:
        /**
         * @brief Default constructor for ROOT I/O
         */
        RayleighScatter() = default;

        /**
         * @brief Construct a Rayleigh scatter object
         * @param position Global position of the Rayleigh scatter
         * @param incoming_momentum Incoming photon momentum vector
         * @param outgoing_momentum Outgoing photon momentum vector
         * @param global_time Global time of the scatter
         * @param track_id Geant4 track ID of the photon
         * @param parent_id Geant4 parent track ID of the photon
         */
        RayleighScatter(ROOT::Math::XYZPoint position,
                        ROOT::Math::XYZVector incoming_momentum,
                        ROOT::Math::XYZVector outgoing_momentum,
                        double global_time,
                        int track_id,
                        int parent_id);

        /**
         * @brief Get the global position of the Rayleigh scatter
         * @return Scatter position
         */
        ROOT::Math::XYZPoint getPosition() const;

        /**
         * @brief Get the incoming photon momentum vector
         * @return Incoming momentum
         */
        ROOT::Math::XYZVector getIncomingMomentum() const;

        /**
         * @brief Get the outgoing photon momentum vector
         * @return Outgoing momentum
         */
        ROOT::Math::XYZVector getOutgoingMomentum() const;

        /**
         * @brief Get the global time of the Rayleigh scatter
         * @return Global time
         */
        double getGlobalTime() const;

        /**
         * @brief Get the Geant4 track ID of the photon
         * @return Geant4 track ID
         */
        int getTrackID() const;

        /**
         * @brief Get the Geant4 parent track ID of the photon
         * @return Geant4 parent track ID
         */
        int getParentID() const;

        /**
         * @brief Print an ASCII representation of RayleighScatter to the given stream
         * @param out Stream to print to
         */
        void print(std::ostream& out) const override;

        void loadHistory() override;
        void petrifyHistory() override;

        /**
         * @brief ROOT class definition
         */
        ClassDefOverride(RayleighScatter, 1); // NOLINT

    private:
        ROOT::Math::XYZPoint position_;
        ROOT::Math::XYZVector incoming_momentum_;
        ROOT::Math::XYZVector outgoing_momentum_;

        double global_time_{};

        int track_id_{};
        int parent_id_{};
    };

    /**
     * @brief Typedef for message carrying Rayleigh scatter objects
     */
    using RayleighScatterMessage = Message<RayleighScatter>;

} // namespace allpix

#endif /* ALLPIX_RAYLEIGH_SCATTER_H */
