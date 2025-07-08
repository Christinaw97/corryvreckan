/**
 * @file
 * @brief File including all current objects
 * @copyright Copyright (c) 2019-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "Cluster.hpp"
#include "MCParticle.hpp"
#include "Multiplet.hpp"
#include "Pixel.hpp"
#include "TimerSignal.hpp"
#include "Track.hpp"
#include "Waveform.hpp"

namespace corryvreckan {
    /**
     * @brief Tuple containing all objects
     */
    using OBJECTS = std::tuple<Cluster, MCParticle, Pixel, TimerSignal, StraightLineTrack, GblTrack, Multiplet, Waveform>;
} // namespace corryvreckan
