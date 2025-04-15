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

#include <TRef.h>
#include <string>
#include <typeindex>

#include <iostream>

#include "Object.hpp"

namespace corryvreckan {
    enum class TimerType {
        NONE,           ///< Indeterminate timer signal
        TRIGGER,        ///< Timer refers to a trigger signal
        POWER_ON,       ///< Timer refers to a powering-on signal
        POWER_OFF,      ///< Timer refers to a power-off signal
        SHUTTER_OPEN,   ///< Timer refers to a shutter opening signal
        SHUTTER_CLOSED, ///< Timer refers to a shutter closing signal
    };

    /**
     * @ingroup Objects
     * @brief Timing signal recorded by a readout system, such as e.g. a trigger
     */
    class TimerSignal : public Object {

    public:
        // Constructors and destructors
        /**
         * @brief Required default constructor
         */
        TimerSignal() = default;

        /**
         * @brief Construct a timer signal without type
         *
         * @param detectorID Name of the detector providing the timer signal
         * @param timestamp Absolute timestamp of the timer signal
         */

        TimerSignal(std::string detectorID, double timestamp) : Object(std::move(detectorID), timestamp){};

        /**
         * @brief Construct timer signal with type
         *
         * @param detectorID Name of the detector providing the timer signal
         * @param timestamp Absolute timestamp of the timer signal
         * @param type Type of the timer signal
         */
        TimerSignal(std::string detectorID, double timestamp, TimerType type)
            : Object(std::move(detectorID), timestamp), type_(type){};

        /**
         * @brief Static member function to obtain base class for storage on the clipboard.
         * This method is used to store objects from derived classes under the typeid of their base classes
         *
         * @warning This function should not be implemented for derived object classes
         *
         * @return Class type of the base object
         */
        static std::type_index getBaseType() { return typeid(TimerSignal); }

        /**
         * @brief Set tag
         *
         * @param tag Tag to set
         */
        void setTag(std::string tag) { tag_ = std::move(tag); }

        /**
         * qbrief Set trigger ID this timer signal should be associated with
         *
         * @param trigger_id Trigger ID to be stored
         */
        void setTriggerID(uint32_t trigger_id) { trigger_id_ = trigger_id; };

        /**
         * @brief Obtain timer signal type
         * @return Timer type
         */
        TimerType getType() const { return type_; }

        /**
         * @brief Obtain timer tag
         * @return Tag of the timer signal
         */
        std::string getTag() const { return tag_; }

        /**
         * @brief Obtain the trigger ID this timer signal is associated to
         * @return Trigger ID
         */
        size_t getTriggerID() const { return trigger_id_; }

        // ROOT I/O class definition - update version number when you change this class!
        ClassDefOverride(TimerSignal, 6);

        void loadHistory() override{};
        void petrifyHistory() override{};

        /**
         * @brief Print an ASCII representation of TimerSignal to the given stream
         * @param out Stream to print to
         */
        void print(std::ostream& out) const override;

    protected:
        TimerType type_{TimerType::NONE};
        std::string tag_{};
        uint32_t trigger_id_{0};
    };

    // Vector type declaration
    using TimerSignalVector = std::vector<std::shared_ptr<TimerSignal>>;
} // namespace corryvreckan

#endif // CORRYVRECKAN_TIMERSIGNAL_H
