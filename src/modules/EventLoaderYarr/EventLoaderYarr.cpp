/**
 * @file
 * @brief Implementation of module EventLoaderYarr
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EventLoaderYarr.h"

using namespace corryvreckan;

EventLoaderYarr::EventLoaderYarr(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), m_detector(detector) {

    // Get the input directory from the configuration
    m_inputDirectory = config_.getPath("input_directory");
    LOG(INFO) << "Detector name: " << m_detector->getName();

}

void EventLoaderYarr::initialize() {
    // Initialise event counter
    m_eventNumber = 0;

    // Locate and open the raw data file
    m_filename = findRawFile(m_inputDirectory, m_detector->getName());
    LOG(INFO) << "Opening file " << m_filename;
    fileHandle.open(m_filename, std::istream::in | std::istream::binary);
    if(!fileHandle.good()){ 
        throw(std::invalid_argument("Path " + m_filename + " could not be opened!"));
    } 
    LOG(DEBUG) << "Opened file " << m_filename;

    
    // Create and configure ROOT Histograms:

    // Hit map
    std::string title = m_detector->getName() + " Hit map";
    hHitMap = new TH2F("hitMap", title.c_str(), m_detector->nPixels().X(), -0.5, m_detector->nPixels().X() - 0.5, m_detector->nPixels().Y(), -0.5, m_detector->nPixels().Y() - 0.5);
    gStyle->SetPalette(kGreyScale);

    // Number of hits throughout a day
    title = m_detector->getName() + " Number of Hits vs. time; time [s]; # hits";
    numHitsVsTime = new TH1F("numHitsVsTime", title.c_str(), 2880, 0, 86400);

    // Number of events vs absolute time (time across whole run)
    title = m_detector->getName() + " Number of Events vs. absolute time; time [s]; # events";
    eventsVsAbsoluteTime = new TH1F("eventsVsAbsoluteTime", title.c_str(), 2880, 0, 86400);
    eventsVsAbsoluteTime->GetXaxis()->SetCanExtend(kTRUE);

    // Number of hits vs absolute time (time across whole run)
    title = m_detector->getName() + " Number of Hits vs. absolute time; time [s]; # hits";
    numHitsVsAbsTime = new TH1F("numHitsVsAbsTime", title.c_str(), 2880, 0, 86400);
    numHitsVsAbsTime->GetXaxis()->SetCanExtend(kTRUE);

}

StatusCode EventLoaderYarr::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Create PixelVector to store hits and buffers to store trigger data
    PixelVector pixels;
    std::vector<uint16_t> triggerL1Ids;           
    std::vector<uint32_t> triggerTimes;
    
    // Read the first trigger header to establish the current trigger window.
    std::streampos headerStartPos = fileHandle.tellg();
    auto firstHeader = readHeader();
    uint16_t current_bcid = firstHeader.bcid;
    triggerL1Ids.push_back(firstHeader.l1id);
    triggerTimes.push_back(firstHeader.time);
    if (firstHeader.numHits > 0) {
        readHits(pixels, firstHeader.time, firstHeader.numHits);
    }

    // Loop through rest of trigger window
    while (fileHandle.peek() != EOF) {
        headerStartPos = fileHandle.tellg();
        auto header = readHeader();
        
        // If the BCID has changed, rewind to before the header so the new event can handle it.
        if (header.bcid != current_bcid) {
            fileHandle.seekg(headerStartPos);
            break;
        }

        triggerL1Ids.push_back(header.l1id);
        triggerTimes.push_back(header.time);
        readHits(pixels, header.time, header.numHits);
    }

    // Compute the absolute event time
    m_absolute_event_time = adjustEventTime(triggerTimes);

    // Get or create the event in the clipboard
    std::shared_ptr<Event> event;
    if(!clipboard->isEventDefined()) {
        event = std::make_shared<Event>(m_absolute_event_time, m_absolute_event_time);
        clipboard->putEvent(event);
    } else {
        event = clipboard->getEvent();
    }

    // Add trigger entries from the buffers to the event
    for (const auto& l1id : triggerL1Ids) {
        event->addTrigger(l1id, static_cast<double>(m_absolute_event_time));
    }

    // Add the pixels to the clipboard
    if(!pixels.empty()) {
        clipboard->putData(pixels, m_detector->getName());
        LOG(DEBUG) << "Added " << pixels.size() << " pixels.";
    }

    // Fill the plots
    eventsVsAbsoluteTime->Fill(static_cast<double>(m_absolute_event_time) / 1000);          // Convert from ms to s
    for(auto& pixel : pixels) {
        hHitMap->Fill(pixel->column(), pixel->row());                                       // Fill the hit map
        numHitsVsTime->Fill(pixel->timestamp() / 1.0e9);                                    // Convert from ns to s
        numHitsVsAbsTime->Fill(static_cast<double>(m_absolute_event_time) / 1000);          // Convert from ms to s
    }
    
    // Finish up
    m_eventNumber++;
    if(fileHandle.peek() == EOF) {
        LOG(STATUS) << "Reached end-of-file. Closing file.";
        fileHandle.close();
        return StatusCode::EndRun;
    }
    return StatusCode::Success;

}

void EventLoaderYarr::finalize(const std::shared_ptr<ReadonlyClipboard>&) { 
    LOG(INFO) << "Analysed " << m_eventNumber << " events"; 
}

// Helper functions:

EventLoaderYarr::TriggerHeader EventLoaderYarr::readHeader() {
    TriggerHeader header;

    fileHandle.read(reinterpret_cast<char*> (&header.tag), sizeof(uint32_t));
    fileHandle.read(reinterpret_cast<char*> (&header.l1id), sizeof(uint16_t));
    fileHandle.read(reinterpret_cast<char*> (&header.bcid), sizeof(uint16_t));
    fileHandle.read(reinterpret_cast<char*> (&header.numHits), sizeof(uint16_t));

    // Derive additional header information from tag:
    header.basetag = static_cast<uint8_t>((header.tag & 252) >> 2);
    header.exttag  = static_cast<uint8_t>(header.tag & 3);
    header.time    = ((header.tag >> 8) & 0xFFFFFF) << 3;

    return header;
}


void EventLoaderYarr::readHits(PixelVector& pixels, uint32_t eventTime, uint16_t nHits) {
    for(uint16_t i = 0; i < nHits; i++) {
        uint16_t col = 0, row = 0, tot = 0;

        fileHandle.read(reinterpret_cast<char*>(&col), sizeof(uint16_t));
        fileHandle.read(reinterpret_cast<char*>(&row), sizeof(uint16_t));
        fileHandle.read(reinterpret_cast<char*>(&tot), sizeof(uint16_t));

        auto pixel = std::make_shared<Pixel>(
            m_detector->getName(),                           // Detector name.
            col,                                             // Column.
            row,                                             // Row.
            tot,                                             // Raw 
            static_cast<double>(tot),                        // Charge: ToT is the charge.
            static_cast<double>(eventTime)*1e6               // Timestamp: every pixel in this event has the same time.
        );

        pixels.push_back(pixel);
    }
}

std::string EventLoaderYarr::findRawFile(const std::string& directory, const std::string& detectorName) {
    // Open the input directory
    DIR* dir = opendir(directory.c_str());

    // Check if the directory exists
    if(dir == nullptr) {
        throw ModuleError("Directory " + directory + " does not exist");
    } else {
        LOG(INFO) << "Found directory " << directory;
    }

    std::vector<std::string> files;
    dirent* entry;

    // Read all .raw files in the directory
    while((entry = readdir(dir))) {
        std::string filename = entry->d_name;
        if(filename.find(".raw") != std::string::npos && filename.find(detectorName) != std::string::npos) {
            files.push_back(directory + "/" + filename);
            LOG(INFO) << "Found YARR data file: " << filename;
        }
    }
    closedir(dir);

    // Check if any files were found
    if(files.empty()) {
        throw ModuleError("No raw data file found for detector " + detectorName + " in " + directory);
    }
    else if (files.size() > 1) {
        throw ModuleError("Multiple raw data files found for detector " + detectorName + " in " + directory);
    }

    return files[0];
}

uint64_t EventLoaderYarr::adjustEventTime(const std::vector<uint32_t>& triggerTimes) {
    // Use the first trigger's time as the current event time.
    uint64_t current_time = (!triggerTimes.empty()) ? triggerTimes.front() : 0;

    if (m_previous_time != 0 && current_time < m_previous_time) {
        m_day_offset++; // increment day rollover counter
        LOG(DEBUG) << "Day rollover detected. Total day offset: " << m_day_offset;
    }

    m_previous_time = current_time;
    m_absolute_event_time = current_time + m_day_offset * 86400000; // 86400000 ms in a day

    return m_absolute_event_time;
}