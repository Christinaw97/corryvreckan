/**
 * @file
 * @brief Implementation of module EventFilter
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "EventFilter.h"

using namespace corryvreckan;

EventFilter::EventFilter(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)) {}

void EventFilter::initialize() {

    config_.setDefault<unsigned>("minTracks", 0);
    config_.setDefault<unsigned>("maxTracks", 100);
    config_.setDefault<unsigned>("minClusters_per_plane", 0);
    config_.setDefault<unsigned>("maxClusters_per_plane", 100);

    minNumberTracks_ = config_.get<unsigned>("minTracks");
    maxNumberTracks_ = config_.get<unsigned>("maxTracks");
    minClustersPerReference_ = config_.get<unsigned>("minClusters_per_plane");
    maxClustersPerReference_ = config_.get<unsigned>("maxClusters_per_plane");
}

StatusCode EventFilter::run(const std::shared_ptr<Clipboard>& clipboard) {

    eventsTotal_++;
    auto numTracks = clipboard->getData<Track>().size();
    if(numTracks > maxNumberTracks_) {
        eventsSkipped_++;
        LOG(TRACE) << "Number of tracks above maximum";
        return StatusCode::DeadTime;
    } else if(numTracks < minNumberTracks_) {
        eventsSkipped_++;
        LOG(TRACE) << "Number of tracks below minimum";
        return StatusCode::DeadTime;
    }

    // Loop over all reference detectors
    for(auto& detector : get_detectors()) {
        // skip duts and auxilliary
        if(detector->isAuxiliary() || detector->isDUT()) {
            continue;
        }
        std::string det = detector->getName();
        // Check if number of Clusters on plane is within acceptance
        auto numClusters = clipboard->getData<Cluster>(det).size();
        if(numClusters > maxClustersPerReference_) {
            eventsSkipped_++;
            LOG(TRACE) << "Number of Clusters on above maximum";
            return StatusCode::DeadTime;
        }
        if(numClusters < minClustersPerReference_) {
            eventsSkipped_++;
            LOG(TRACE) << "Number of Clusters on below minimum";
            return StatusCode::DeadTime;
        }
    }
    return StatusCode::Success;
}

void EventFilter::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    LOG(STATUS) << "Skipped " << eventsSkipped_ << " events of " << eventsTotal_;
}
