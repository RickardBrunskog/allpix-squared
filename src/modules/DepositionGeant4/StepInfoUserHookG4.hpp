/**
 * @file
 * @brief Defines a user hook for Geant4 stepping action
 *
 * @copyright Copyright (c) 2023-2025 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef StepInfoUserHookG4_H
#define StepInfoUserHookG4_H 1

#include <cstddef>
#include <vector>

#include "G4Step.hh"
#include "G4UserSteppingAction.hh"

#include "objects/RayleighScatter.hpp"

namespace allpix {

    class Event;
    class Messenger;
    class Module;

    /**
     * @brief Allows access to the info of each Geant4 step
     */
    class StepInfoUserHookG4 : public G4UserSteppingAction {
    public:
        /**
         * @brief Default constructor
         */
        explicit StepInfoUserHookG4() = default;

        /**
         * @brief Default destructor
         */
        ~StepInfoUserHookG4() override = default;

        /**
         * @brief Called for every step in Geant4
         * @param aStep The pointer to the G4Step for which this routine is called
         */
        void UserSteppingAction(const G4Step* aStep) override;

        /**
         * @brief Clear stored Rayleigh scatter information for the current event
         */
        static void clearRayleighScatterInfo();

        /**
         * @brief Get the number of stored Rayleigh scatters for the current event
         * @return Number of stored Rayleigh scatters
         */
        static std::size_t getStoredRayleighScatters();

        /**
         * @brief Dispatch stored Rayleigh scatters as an Allpix message
         * @param module Module dispatching the message
         * @param messenger Messenger to dispatch to
         * @param event Current Allpix event
         */
        static void dispatchRayleighScatterMessage(Module* module, Messenger* messenger, Event* event);

    private:
        static thread_local std::vector<RayleighScatter> rayleigh_scatters_;
    };

} // namespace allpix

#endif /* StepInfoUserHookG4_H */