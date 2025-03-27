/**
 * @file
 * @brief Implementation of module EventLoaderHDF5
 *
 * @copyright Copyright (c) 2023-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EventLoaderHDF5.h"

namespace corryvreckan {

    EventLoaderHDF5::EventLoaderHDF5(Configuration& config, std::shared_ptr<Detector> detector)
        : Module(config, detector), m_detector(detector) {
        h5_datatype.insertMember("column", HOFFSET(Hit, column), H5::PredType::STD_U16LE);
        h5_datatype.insertMember("row", HOFFSET(Hit, row), H5::PredType::STD_U16LE);
        h5_datatype.insertMember("charge", HOFFSET(Hit, charge), H5::PredType::STD_U8LE);
        h5_datatype.insertMember("timestamp", HOFFSET(Hit, timestamp), H5::PredType::STD_U64LE);
        h5_datatype.insertMember("trigger_number", HOFFSET(Hit, trigger_number), H5::PredType::STD_U32LE);

        m_fileName = config.getPath("filename");
        m_datasetName = config.get<std::string>("dataset_name", "Hits");
        m_bufferDepth = config.get<hsize_t>("buffer_depth", 100000);
        m_sync_by_trigger = config.get<bool>("sync_by_trigger", false);
        m_eventLength = config.get<double>("event_length", Units::get<double>(1.0, "us"));
        m_timestampShift = config.get<double>("timestamp_shift", 0);
        m_triggerShift = config.get<uint32_t>("trigger_shift", 0);
    }

    void EventLoaderHDF5::initialize() {
        // Open the file
        try {
            m_file = H5::H5File(m_fileName, H5F_ACC_RDONLY);
            m_dataset = m_file.openDataSet(m_datasetName);
        } catch(const std::exception& ex) {
            LOG(ERROR) << "Failed to open file.";
            throw;
        }

        if(!m_detector->isAuxiliary()) {
            // Initialize hitmap and charge histograms
            hHitMap = new TH2F("hitMap",
                               "Hit Map",
                               m_detector->nPixels().X(),
                               -0.5,
                               m_detector->nPixels().X() - 0.5,
                               m_detector->nPixels().Y(),
                               -0.5,
                               m_detector->nPixels().Y() - 0.5);
            hTotMap = new TH2F("totMap",
                               "ToT Map",
                               m_detector->nPixels().X(),
                               -0.5,
                               m_detector->nPixels().X() - 0.5,
                               m_detector->nPixels().Y(),
                               -0.5,
                               m_detector->nPixels().Y() - 0.5);
            hPixelToT = new TH1F("pixelToT", "Pixel ToT", 200, -0.5, 199.5);
        }

        // TODO: only define those if event is not defined yet. How to find out here?
        std::string title =
            "Corryvreckan event start times (placed on clipboard); Corryvreckan event start time [ms];# entries";
        hClipboardEventStart = new TH1D("clipboardEventStart", title.c_str(), 3e6, 0, 3e3);

        title = "Corryvreckan event start times (placed on clipboard); Corryvreckan event start time [s];# entries";
        hClipboardEventStart_long = new TH1D("clipboardEventStart_long", title.c_str(), 3e6, 0, 3e3);

        title = "Corryvreckan event end times (placed on clipboard); Corryvreckan event end time [ms];# entries";
        hClipboardEventEnd = new TH1D("clipboardEventEnd", title.c_str(), 3e6, 0, 3e3);

        title = "Corryvreckan event durations (on clipboard); Corryvreckan event duration [ms];# entries";
        hClipboardEventDuration = new TH1D("clipboardEventDuration", title.c_str(), 3e6, 0, 3e3);

        m_start_record = 0;

        f_dataspace = m_dataset.getSpace();
        f_total_records = static_cast<hsize_t>(f_dataspace.getSimpleExtentNpoints());
        LOG(DEBUG) << "Total number of records " << f_total_records;
    }

    StatusCode EventLoaderHDF5::run(const std::shared_ptr<Clipboard>& clipboard) {
        PixelVector deviceData;

        // Load data directly into the vector
        bool data = loadData(clipboard, deviceData);

        if(data) {
            LOG(DEBUG) << "Loaded " << deviceData.size() << " pixels for device " << m_detector->getName();
            clipboard->putData(deviceData, m_detector->getName());
        }

        LOG(DEBUG) << clipboard->countObjects<Pixel>() << " objects on the clipboard";
        if(m_buffer.empty() && (m_start_record == f_total_records)) {
            return StatusCode::EndRun;
        }

        return StatusCode::Success;
    }

    bool EventLoaderHDF5::loadData(const std::shared_ptr<Clipboard>& clipboard, PixelVector& deviceData_) {
        // Ensure that m_buffer is filled with data
        fillBuffer();

        std::string detectorID = m_detector->getName();

        while(!m_buffer.empty()) {
            auto hit = m_buffer.top();

            double shiftedTimestamp = static_cast<double>(hit->timestamp) + m_timestampShift;
            uint32_t shiftedTriggerId = hit->trigger_number + m_triggerShift;

            // Check if an event is defined or if we need to create it:
            if(!clipboard->isEventDefined()) {
                double event_start = shiftedTimestamp;
                double event_end = event_start + m_eventLength;
                LOG(DEBUG) << "Defining Corryvreckan event: " << Units::display(event_start, {"us", "ns"}) << " - "
                           << Units::display(event_end, {"us", "ns"}) << ", length "
                           << Units::display(event_end - event_start, {"us", "ns"});
                clipboard->putEvent(std::make_shared<Event>(event_start, event_end));
                clipboard->getEvent()->addTrigger(shiftedTriggerId,
                                                  event_start); // TODO: decide where to put trigger inside the event? Maybe
                                                                // also in case of already defined events?
                hClipboardEventStart->Fill(static_cast<double>(Units::convert(event_start, "ms")));
                hClipboardEventStart_long->Fill(static_cast<double>(Units::convert(event_start, "s")));
                hClipboardEventEnd->Fill(static_cast<double>(Units::convert(event_end, "ms")));
                hClipboardEventDuration->Fill(static_cast<double>(
                    Units::convert(clipboard->getEvent()->end() - clipboard->getEvent()->start(), "ms")));
            } else {
                LOG(DEBUG) << "Corryvreckan event found on clipboard: "
                           << Units::display(clipboard->getEvent()->start(), {"us", "ns"}) << " - "
                           << Units::display(clipboard->getEvent()->end(), {"us", "ns"})
                           << ", length: " << Units::display(clipboard->getEvent()->duration(), {"us", "ns"});
            }

            auto event = clipboard->getEvent();
            Event::Position position = getPosition(event, hit);

            if(position == Event::Position::AFTER) {
                LOG(DEBUG) << "Stopping processing event, pixel is after event window ("
                           << Units::display(shiftedTimestamp, {"s", "us", "ns"}) << " > "
                           << Units::display(event->end(), {"s", "us", "ns"}) << ")";
                break;
            } else if(position == Event::Position::BEFORE) {
                LOG(TRACE) << "Skipping pixel, is before event window ("
                           << Units::display(shiftedTimestamp, {"s", "us", "ns"}) << " < "
                           << Units::display(event->start(), {"s", "us", "ns"}) << ")";
                m_buffer.pop();
            } else {
                LOG(DEBUG) << "Position is DURING";
                if(!m_detector->isAuxiliary()) {
                    double pixel_timestamp;
                    LOG(DEBUG) << "Loaded pixel (" << hit->column << ", " << hit->row << ")";
                    if(m_sync_by_trigger) {
                        pixel_timestamp =
                            event->getTriggerTime(shiftedTriggerId) + m_timestampShift; // Use trigger time as pixel time
                    } else {
                        pixel_timestamp = shiftedTimestamp;
                    }
                    auto pixel = std::make_shared<Pixel>(
                        m_detector->getName(), hit->column, hit->row, hit->charge, hit->charge, pixel_timestamp);
                    deviceData_.push_back(pixel);
                    hHitMap->Fill(pixel->column(), pixel->row());
                    hTotMap->Fill(pixel->column(), pixel->row(), pixel->raw());
                    hPixelToT->Fill(pixel->raw());
                }
                m_buffer.pop();
            }
            // Refill buffer for next iteration
            fillBuffer();
        }

        if(deviceData_.empty()) {
            return false;
        }

        return true;
    }

    std::vector<EventLoaderHDF5::Hit> EventLoaderHDF5::readChunk() {
        hsize_t num_records_to_read = std::min(m_bufferDepth, f_total_records - m_start_record);
        H5::DataSpace mem_space = H5::DataSpace(1, &num_records_to_read);
        std::vector<EventLoaderHDF5::Hit> chunk(num_records_to_read);

        // Select memory space within the file to read
        f_dataspace.selectHyperslab(H5S_SELECT_SET, &num_records_to_read, &m_start_record, nullptr, nullptr);
        m_dataset.read(chunk.data(), h5_datatype, mem_space, f_dataspace);
        m_start_record += num_records_to_read;

        return chunk;
    }

    void EventLoaderHDF5::fillBuffer() {
        // Fill buffer only if it is empty and there are records left in the file
        if(m_buffer.empty() && (m_start_record != f_total_records)) {
            std::vector<Hit> chunk = readChunk();

            // Add the elements of chunk to the buffer
            for(auto hit : chunk) {
                m_buffer.push(std::make_shared<EventLoaderHDF5::Hit>(hit));
            }
        }
    }

    Event::Position EventLoaderHDF5::getPosition(const std::shared_ptr<Event>& event,
                                                 const std::shared_ptr<Hit>& hit) const {

        uint32_t shiftedTriggerId = hit->trigger_number + m_triggerShift;
        double shiftedTimestamp = static_cast<double>(hit->timestamp) + m_timestampShift;

        if(m_sync_by_trigger) {
            const auto trigger_position = event->getTriggerPosition(shiftedTriggerId);
            LOG(DEBUG) << "Corryvreckan event with trigger id " << shiftedTriggerId << " has trigger time at "
                       << Units::display(event->getTriggerTime(shiftedTriggerId), {"s", "us", "ns"});
            if(trigger_position == Event::Position::BEFORE) {
                LOG(DEBUG) << "(Shifted) trigger ID " << shiftedTriggerId
                           << " is before triggers registered in Corryvreckan event";
                // LOG(DEBUG) << "(Shifted) Trigger ID " << trigger_after_shift
            } else if(trigger_position == Event::Position::AFTER) {
                LOG(DEBUG) << "(Shifted) trigger ID " << shiftedTriggerId
                           << " is after triggers registered in Corryvreckan event";
            } else if(trigger_position == Event::Position::UNKNOWN) {
                LOG(DEBUG) << "(Shifted) trigger ID " << shiftedTriggerId
                           << " is within Corryvreckan event range but not registered";
            }
            return trigger_position;
        } else {
            return event->getTimestampPosition(shiftedTimestamp);
        }
    }

} // namespace corryvreckan
