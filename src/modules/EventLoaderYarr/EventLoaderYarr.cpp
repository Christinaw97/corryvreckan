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
 
 constexpr double MS_TO_NS = 1e6;
 constexpr double NS_TO_S = 1e-9;
 constexpr double DAY_IN_NS = 86400000000000;
 
 EventLoaderYarr::EventLoaderYarr(Configuration& config, std::shared_ptr<Detector> detector)
     : Module(config, detector), detector_(detector) {
     // Get config parameters
     input_directory_ = config_.getPath("input_directory");
     trigger_tag_timing_ = config_.get<bool>("trigger_tag_timing", false);
 
     if(trigger_tag_timing_)
         LOG(INFO) << "Using tag timestamp for event time.";
 }
 
 void EventLoaderYarr::initialize() {
     event_number_ = 0;
     file_name_ = find_raw_file(input_directory_, detector_->getName());
     file_handle_.open(file_name_, std::istream::in | std::istream::binary);
     if(!file_handle_.good())
         throw ModuleError("Cannot open file: " + file_name_);
 
     // Create and configure ROOT Histograms:
 
     // Hit map
     std::string title = detector_->getName() + " Hit map";
     hHitMap = new TH2F("hitMap", (detector_->getName() + " Hit map").c_str(), detector_->nPixels().X(), -0.5,
                        detector_->nPixels().X() - 0.5, detector_->nPixels().Y(), -0.5, detector_->nPixels().Y() - 0.5);
 
     if(trigger_tag_timing_) {
         // Number of events vs time
         title = detector_->getName() + " Number of Events vs. Tag Timestamp; time [s]; # events";
         hEventsVsTagTime = new TH1F("eventsVsTimestamp", title.c_str(), 2880, 0, 86400);
         hEventsVsTagTime->GetXaxis()->SetCanExtend(kTRUE);
 
         // Number of hits vs time
         title = detector_->getName() + " Number of Hits vs. Tag Timestamp; time [s]; # hits";
         hNumHitsVsTagTime = new TH1F("numHitsVsTimestamp", title.c_str(), 2880, 0, 86400);
         hNumHitsVsTagTime->GetXaxis()->SetCanExtend(kTRUE);
     }
 }
 
 StatusCode EventLoaderYarr::run(const std::shared_ptr<Clipboard>& clipboard) {
     PixelVector pixels;
     std::vector<uint16_t> triggerL1Ids;
     std::vector<double> triggerTimes;
     std::streampos headerStartPos = file_handle_.tellg();
 
     auto firstHeader = read_header();
     double first_timestamp = (trigger_tag_timing_ ? adjust_time(firstHeader.time) : firstHeader.time);
     triggerL1Ids.push_back(firstHeader.l1id);
     triggerTimes.push_back(first_timestamp);
     if(firstHeader.bcid != event_number_)
         LOG(WARNING) << "BCID vs Event Number Desynchronization: " << firstHeader.bcid << " vs. " << event_number_;
     if(firstHeader.numHits > 0)
         read_hits(pixels, first_timestamp, firstHeader.numHits);
 
     // Process subsequent headers within the same BCID window
     double lastTimestamp = first_timestamp;
     while(file_handle_.peek() != EOF) {
         headerStartPos = file_handle_.tellg();
         auto header = read_header();
         if(header.bcid != firstHeader.bcid) {
             file_handle_.seekg(headerStartPos);
             break;
         }
 
         double event_time_ns = (trigger_tag_timing_ ? (first_timestamp + (header.l1id * 25)) : header.time);
         triggerL1Ids.push_back(header.l1id);
         triggerTimes.push_back(event_time_ns);
         read_hits(pixels, event_time_ns, header.numHits);
         lastTimestamp = event_time_ns;
     }
 
     // Get or create the event in the clipboard
     std::shared_ptr<Event> event;
     if(!clipboard->isEventDefined()) {
         event = std::make_shared<Event>(first_timestamp, lastTimestamp);
         clipboard->putEvent(event);
     } else {
         event = clipboard->getEvent();
         for(const auto& triggerTime : triggerTimes) {
             const auto position = event->getTimestampPosition(triggerTime);
             if(position != Event::Position::DURING) {
                 LOG(WARNING) << "Event timestamp (" << first_timestamp << ") is not in the expected position between "
                              << event->start() << " and " << event->end() << ".";
             }
         }
     }
 
     // Add trigger entries and pixels to the event
     for(size_t i = 0; i < triggerL1Ids.size(); i++)
         event->addTrigger(triggerL1Ids[i], triggerTimes[i]);
     if(!pixels.empty()) {
         clipboard->putData(pixels, detector_->getName());
         LOG(DEBUG) << "Added " << pixels.size() << " pixels to event " << firstHeader.bcid;
     }
 
     // Fill ROOT histograms
     if(trigger_tag_timing_)
         hEventsVsTagTime->Fill(first_timestamp * NS_TO_S);
     for(const auto& pixel : pixels) {
         hHitMap->Fill(pixel->column(), pixel->row());
         if(trigger_tag_timing_)
             hNumHitsVsTagTime->Fill(pixel->timestamp() * NS_TO_S);
     }
 
     // Finish up
     event_number_++;
     if(file_handle_.eof()) {
         LOG(STATUS) << "Reached end-of-file. Closing file.";
         file_handle_.close();
         return StatusCode::EndRun;
     }
     return StatusCode::Success;
 }
 
 void EventLoaderYarr::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
     LOG(INFO) << "Analysed " << event_number_ << " events";
 }
 
 // Helper functions:
 
 EventLoaderYarr::TriggerHeader EventLoaderYarr::read_header() {
     TriggerHeader header;
 
     file_handle_.read(reinterpret_cast<char*>(&header.tag), sizeof(uint32_t));
     file_handle_.read(reinterpret_cast<char*>(&header.l1id), sizeof(uint16_t));
     file_handle_.read(reinterpret_cast<char*>(&header.bcid), sizeof(uint16_t));
     file_handle_.read(reinterpret_cast<char*>(&header.numHits), sizeof(uint16_t));
 
     // Derive additional header information:
     if(trigger_tag_timing_) {
         header.time = (static_cast<double>(((header.tag >> 8) & 0xFFFFFF) << 3)) * MS_TO_NS;
     } else {
         header.time = (header.bcid * 0.025 * MS_TO_NS) + (header.l1id * 25); // 25 ns per L1ID
     }
 
     return header;
 }
 
 void EventLoaderYarr::read_hits(PixelVector& pixels, double eventTime, uint16_t nHits) {
     for(uint16_t i = 0; i < nHits; i++) {
         uint16_t col = 0, row = 0, tot = 0;
 
         file_handle_.read(reinterpret_cast<char*>(&col), sizeof(uint16_t));
         file_handle_.read(reinterpret_cast<char*>(&row), sizeof(uint16_t));
         file_handle_.read(reinterpret_cast<char*>(&tot), sizeof(uint16_t));
 
         pixels.emplace_back(std::make_shared<Pixel>(detector_->getName(),     // Detector name.
                                                     col,                      // Column.
                                                     row,                      // Row.
                                                     tot,                      // Raw
                                                     static_cast<double>(tot), // Charge: ToT is the charge.
                                                     eventTime // Timestamp: every pixel in this trigger has the same time.
                                                     ));
     }
 }
 
 std::string EventLoaderYarr::find_raw_file(const std::string& directory, const std::string& detectorName) {
     DIR* dir = opendir(directory.c_str());
     if(dir == nullptr)
         throw ModuleError("Directory " + directory + " does not exist");
     std::vector<std::string> files;
     dirent* entry;
 
     // Read all .raw files in the directory
     while((entry = readdir(dir))) {
         std::string filename = entry->d_name;
         if(filename.find(".raw") != std::string::npos && filename.find(detectorName) != std::string::npos) {
             files.push_back(directory + "/" + filename);
             LOG(INFO) << "Found a data file named " << filename << " for detector " << detectorName;
         }
     }
     closedir(dir);
     if(files.empty()) {
         throw ModuleError("No raw data file found for detector " + detectorName + " in " + directory);
     } else if(files.size() > 1) {
         throw ModuleError("Multiple raw data files found for detector " + detectorName + " in " + directory);
     }
 
     return files[0];
 }
 
 double EventLoaderYarr::adjust_time(const double& current_time) {
     if(current_time < previous_time_) {
         day_offset_++;
         LOG(DEBUG) << "Day rollover detected. Total day offset: " << day_offset_;
     }
     previous_time_ = current_time;
 
     return current_time + day_offset_ * DAY_IN_NS;
 }