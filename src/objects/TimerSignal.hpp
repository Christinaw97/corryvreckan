/**
 * @file
 * @brief Definition of SPIDR signal object
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRYVRECKAN_TIMERSIGNAL_H
#define CORRYVRECKAN_TIMERSIGNAL_H 1

#include <string>
#include <typeindex>

#include "Object.hpp"

namespace corryvreckan {
    /**
     * @ingroup Objects
     * @brief Timing signal recorded by a readout system, such as e.g. a trigger
     */
    class TimerSignal : public Object {

    public:
        // Constructors and destructors
        TimerSignal() {};
        TimerSignal(double timestamp) : Object(timestamp) {};

        void setTag(std::string tag) { tag_ = std::move(tag); }

        void setTriggerID(uint32_t trigger_id) { trigger_id_ = trigger_id; };

        /**
         * @brief Static member function to obtain base class for storage on the clipboard.
         * This method is used to store objects from derived classes under the typeid of their base classes
         *
         * @warning This function should not be implemented for derived object classes
         *
         * @return Class type of the base object
         */
        static std::type_index getBaseType() { return typeid(TimerSignal); }

        std::string getTag() const { return tag_; }
        size_t getTriggerID() const { return trigger_id_; }

        void loadHistory() override {};
        void petrifyHistory() override {};

        // ROOT I/O class definition - update version number when you change this class!
        ClassDefOverride(TimerSignal, 4);

    protected:
        std::string tag_;
        uint32_t trigger_id_;
    };

    // Vector type declaration
    using TimerSignalVector = std::vector<std::shared_ptr<TimerSignal>>;
} // namespace corryvreckan

#endif // CORRYVRECKAN_TIMERSIGNAL_H
