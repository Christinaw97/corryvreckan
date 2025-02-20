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
    m_triggerMultiplicity = config_.get<int>("trigger_multiplicity", 17);

    LOG(INFO) << "Detector name: " << m_detector->getName();
    LOG(INFO) << "Trigger multiplicity has been set to " << m_triggerMultiplicity;

}

void EventLoaderYarr::initialize() {

    // Initialise member variables
    m_eventNumber = 0;

    // Open the input directory
    DIR* directory = opendir(m_inputDirectory.c_str());
    if(directory == nullptr) {
        throw ModuleError("Directory " + m_inputDirectory + " does not exist");
    } else {
        LOG(INFO) << "Found directory " << m_inputDirectory;
    }


    // Read all .raw files in the directory
    dirent* entry;

    // Buffer for file names:
    std::vector<std::string> m_files;

    // Grab the correct data file for the detector
    while((entry = readdir(directory))) {
        std::string filename = entry->d_name;
        if(filename.find(".raw") != std::string::npos && filename.find(m_detector->getName()) != std::string::npos) {
            m_files.push_back(m_inputDirectory + "/" + filename);
            LOG(INFO) << "Found YARR data file: " << filename;
        }
    }
    closedir(directory);
    if(m_files.empty()) {
        throw ModuleError("No raw data file found for detector " + m_detector->getName() + " in input directory.");
    }
    else if (m_files.size() > 1) {
        throw ModuleError("Multiple raw data files found for detector " + m_detector->getName() + " in input directory.");
    }
    filename_ = m_files[0];


    // Open the file
    fileHandle.open(filename_, std::istream::in | std::istream::binary);
    filePos = fileHandle.tellg();
    if(!fileHandle.good()){ 
        throw(std::invalid_argument("Path " + filename_ + " could not be opened!"));
    } else {
        LOG(DEBUG) << "Opened file " << filename_;
    }


    // Create hitmap for detector
    std::string title = m_detector->getName() + " Hit map";
    gStyle->SetPalette(kGreyScale);
    hHitMap = new TH2F("hitMap", title.c_str(), m_detector->nPixels().X(), -0.5, m_detector->nPixels().X() - 0.5, m_detector->nPixels().Y(), -0.5, m_detector->nPixels().Y() - 0.5);
    numHitsVsTime = new TH1F("numHitsVsTime", "Number of Hits vs. time; time [s]; # hits", 2880, 0, 86400);

}

StatusCode EventLoaderYarr::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Create PixelVector to store hits
    PixelVector pixels;

    // Buffers to store trigger data from the trigger window
    std::vector<uint16_t> triggerL1Ids;           
    std::vector<uint32_t> triggerTimes;        

    // Loop through trigger window
    for(int i = 0; i < m_triggerMultiplicity; i++) {
        
        // Read the trigger
        this_tag = 0;
        this_l1id = 0;
        this_bcid = 0;
        this_t_hits = 0;

        // Read the header of the event
        readHeader();
        this_basetag = (this_tag & 252) >> 2;
        this_exttag = (this_tag & 3);
        this_time = ((this_tag >> 8) & (0xFFFFFF)) << 3;

        triggerL1Ids.push_back(this_l1id);
        triggerTimes.push_back(this_time);
        
        // Read the hits of the event
        readHits(pixels);
    }

    m_event_time = (!triggerTimes.empty()) ? triggerTimes.front() : 0;

    std::shared_ptr<Event> event;

    if(!clipboard->isEventDefined()) {
        event = std::make_shared<Event>(m_event_time, m_event_time);
        clipboard->putEvent(event);
    } else {
        event = clipboard->getEvent();
    }

    // Add all trigger entries from the buffers to the event.
    for (size_t i = 0; i < triggerL1Ids.size(); ++i) {
        event->addTrigger(triggerL1Ids[i], triggerTimes[i]);
    }

    if(!pixels.empty()) {
        clipboard->putData(pixels, m_detector->getName());
        LOG(DEBUG) << "Added " << pixels.size() << " pixels.";
    }

    // Fill the plots
    for(auto& pixel : pixels) {
        hHitMap->Fill(pixel->column(), pixel->row());
        numHitsVsTime->Fill(pixel->timestamp() / 1.0e9);
    }

    // Increment event counter
    m_eventNumber++;

    if(fileHandle.peek() == EOF) {
        LOG(INFO) << "Reached end-of-file. Closing file.";
        fileHandle.close();

        // Finish the remaining modules for this event and then end the run
        return StatusCode::EndRun;
    }
    else{
        // Continue with next modules and then next event
        return StatusCode::Success;
    }

}

void EventLoaderYarr::finalize(const std::shared_ptr<ReadonlyClipboard>&) { 

    LOG(STATUS) << "Analysed " << m_eventNumber << " events"; 
    LOG(DEBUG) << "Last BCID: " << this_bcid;
}


void EventLoaderYarr::readHeader() {
    fileHandle.read(reinterpret_cast<char*> (&this_tag), sizeof(uint32_t));
    fileHandle.read(reinterpret_cast<char*> (&this_l1id), sizeof(uint16_t));
    fileHandle.read(reinterpret_cast<char*> (&this_bcid), sizeof(uint16_t));
    fileHandle.read(reinterpret_cast<char*> (&this_t_hits), sizeof(uint16_t));
}


void EventLoaderYarr::readHits(PixelVector& pixels) {
    for(uint16_t i = 0; i < this_t_hits; i++) {
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
            static_cast<double>(this_time)*1000000           // Timestamp: every pixel in this event has the same time.
        );

        pixels.push_back(pixel);
    }
}