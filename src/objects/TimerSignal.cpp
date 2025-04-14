/**
 * @file
 * @brief Implementation of MCParticle object
 *
 * @copyright Copyright (c) 217-2018-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "TimerSignal.hpp"


using namespace corryvreckan;


void TimerSignal::print(std::ostream& out) const {
    out << "TimerSignal " << this->timestamp() << ", " << this->trigger_id_ << ", " << this->tag_;
}
