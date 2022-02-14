/**
 * @file
 * @brief Implementation of module EventFilter
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "FilterEvents.h"

using namespace corryvreckan;

FilterEvents::FilterEvents(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)) {}

void FilterEvents::initialize() {

    config_.setDefault<unsigned>("min_tracks", 0);
    config_.setDefault<unsigned>("max_tracks", 100);
    config_.setDefault<unsigned>("min_clusters_per_plane", 0);
    config_.setDefault<unsigned>("max_clusters_per_plane", 100);
    config_.setDefaultMap<std::string, std::string>("filter_tags", std::map<std::string, std::string>{});

    min_number_tracks_ = config_.get<unsigned>("min_tracks");
    max_number_tracks_ = config_.get<unsigned>("max_tracks");
    min_clusters_per_reference_ = config_.get<unsigned>("min_clusters_per_plane");
    max_clusters_per_reference_ = config_.get<unsigned>("max_clusters_per_plane");
    tag_filters_ = config_.getMap<std::string, std::string>("filter_tags", std::map<std::string, std::string>{});

    hFilter_ = new TH1F("FilteredEvents", "Events filtered;events", 6, 0.5, 6.5);
    hFilter_->GetXaxis()->SetBinLabel(1, "Events");
    hFilter_->GetXaxis()->SetBinLabel(2, "Too few tracks");
    hFilter_->GetXaxis()->SetBinLabel(3, "Too many tracks");
    hFilter_->GetXaxis()->SetBinLabel(4, "Too few clusters");
    hFilter_->GetXaxis()->SetBinLabel(5, "Too many clusters");
    hFilter_->GetXaxis()->SetBinLabel(6, "Events passed ");
}

StatusCode FilterEvents::run(const std::shared_ptr<Clipboard>& clipboard) {

    hFilter_->Fill(1); // number of events
    auto status = filter_tracks(clipboard) ? StatusCode::DeadTime : StatusCode::Success;
    status = filter_cluster(clipboard) ? StatusCode::DeadTime : status;
    status = filter_tags(clipboard) ? StatusCode::DeadTime : status;

    if(status == StatusCode::Success) {
        hFilter_->Fill(6);
    }
    return status;
}

void FilterEvents::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    LOG(STATUS) << "Skipped " << hFilter_->GetBinContent(1) << " events. Events passed " << hFilter_->GetBinContent(6);
}

bool FilterEvents::filter_tracks(const std::shared_ptr<Clipboard>& clipboard) {
    auto num_tracks = clipboard->getData<Track>().size();
    if(num_tracks > max_number_tracks_) {
        hFilter_->Fill(2); // too many tracks
        LOG(TRACE) << "Number of tracks above maximum";
        return true;
    } else if(num_tracks < min_number_tracks_) {
        hFilter_->Fill(3); //  too few tracks
        LOG(TRACE) << "Number of tracks below minimum";
        return true;
    }
    return false;
}

bool FilterEvents::filter_cluster(const std::shared_ptr<Clipboard>& clipboard) {
    // Loop over all reference detectors
    for(auto& detector : get_regular_detectors(false)) {
        std::string det = detector->getName();
        // Check if number of Clusters on plane is within acceptance
        auto num_clusters = clipboard->getData<Cluster>(det).size();
        if(num_clusters > max_clusters_per_reference_) {
            hFilter_->Fill(4); //  too many clusters
            LOG(TRACE) << "Number of Clusters on above maximum";
            return true;
        }
        if(num_clusters < min_clusters_per_reference_) {
            hFilter_->Fill(5); //  too few clusters
            LOG(TRACE) << "Number of Clusters on below minimum";
            return true;
        }
    }
    return false;
}

bool FilterEvents::is_tag_filter_passed(const std::string& tag_value, const std::string& tag_filter) {
    // locate range brackets if they exist
    std::size_t open_bracket_pos = tag_filter.find("[");
    std::size_t close_bracket_pos = tag_filter.find("]");
    if((open_bracket_pos != std::string::npos) && (close_bracket_pos != std::string::npos)) {
        // found range brackets, now fetch range values
        std::vector<double> range_values = corryvreckan::split<double>(
            tag_filter.substr(open_bracket_pos + 1, close_bracket_pos - open_bracket_pos - 1), ":");
        if(range_values.size() > 2) {
            throw std::invalid_argument("invalid key value : tag range should hold two values in brackets, separated by a "
                                        "semicolon. Check for extra semicolon");
        }
        double value = corryvreckan::from_string<double>(tag_value);
        double min_value = range_values.at(0);
        double max_value = range_values.at(1);
        if(value > max_value) {
            LOG(TRACE) << "Tag value above maximum";
            return false;
        }
        if(value < min_value) {
            LOG(TRACE) << "Tag value below minimum";
            return false;
        }
        return true;
    } else {
        std::vector<std::string> tag_filter_values = corryvreckan::split<std::string>(tag_filter, ",");
        for(const auto& filter_value : tag_filter_values) {
            if(tag_value == filter_value) {
                return true;
            }
        }
        LOG(TRACE) << "Tag value different from required";
        return false;
    }
}

bool FilterEvents::filter_tags(const std::shared_ptr<Clipboard>& clipboard) {
    auto event = clipboard->getEvent();
    for(auto& tag_filter_pair : tag_filters_) {
        auto tag_filter_key = tag_filter_pair.first;
        auto tag_filter_value = tag_filter_pair.second;
        try {
            auto tag_value = event->getTag(tag_filter_key);
            LOG(TRACE) << "Applying filter " << tag_filter_value << "  to tag " << tag_filter_key << " with value "
                       << tag_value;
            if(tag_value.empty()) {
                return true;
            } else {
                return !is_tag_filter_passed(event->getTag(tag_filter_key), tag_filter_value);
            }
        } catch(std::out_of_range& e) {
            throw MissingKeyError(tag_filter_key, config_.getName());
        } catch(std::invalid_argument& e) {
            throw InvalidKeyError(tag_filter_key, config_.getName(), tag_filter_value, typeid(std::string), e.what());
        } catch(std::overflow_error& e) {
            throw InvalidKeyError(tag_filter_key, config_.getName(), tag_filter_value, typeid(std::string), e.what());
        }
    }
    return false;
}
