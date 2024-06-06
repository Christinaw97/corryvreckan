/**
 * @file
 * @brief Implementation of module EventDisplay
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EventDisplay.h"

#include <TProfile2D.h>

using namespace corryvreckan;

EventDisplay::EventDisplay(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), detector_(detector) {}

StatusCode EventDisplay::run(const std::shared_ptr<Clipboard>& clipboard) {

    auto pixels = clipboard->getData<Pixel>(detector_->getName());
    if(pixels.empty()) {
        LOG(DEBUG) << "Detector " << detector_->getName() << " does not have any pixels on the clipboard";
        return StatusCode::Success;
    }

    std::string title = "map_event_" + std::to_string(event_number_);
    auto* histogram = new TProfile2D(title.c_str(),
                                     "rawValues; column; row; raw values",
                                     detector_->nPixels().X(),
                                     -0.5,
                                     detector_->nPixels().X() - 0.5,
                                     detector_->nPixels().Y(),
                                     -0.5,
                                     detector_->nPixels().Y() - 0.5);

    for(const auto& pixel : pixels) {
        histogram->Fill(pixel->column(), pixel->row(), pixel->raw());
    }

    histogram->Write();
    event_number_++;

    // Return value telling analysis to keep running
    return StatusCode::Success;
}
