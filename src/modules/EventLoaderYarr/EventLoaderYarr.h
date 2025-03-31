/**
 * @file
 * @brief Definition of module EventLoaderYarr
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

 #include <TCanvas.h>
 #include <TH1F.h>
 #include <TH2F.h>
 #include <thread>
 #include <chrono>
 #include <iostream>
 #include <queue>
 #include <dirent.h>
 #include "core/module/Module.hpp"
 #include "objects/Cluster.hpp"
 #include "objects/Pixel.hpp"
 #include "objects/Track.hpp"
 
 namespace corryvreckan {
     /** @ingroup Modules
      * @brief Module to load ATLAS ITkPixV2 YARR events from a raw file.
      *
      * This module reads a .raw file from a specified input directory (and detector name)
      * and extracts trigger and hit data to be added to the clipboard. It also creates ROOT plots.
      */
     class EventLoaderYarr : public Module {
 
     public:
         /**
          * @brief Constructor: retrieves configuration values (input directory, trigger multiplicity)
          * @param config Configuration object for this module as retrieved from the steering file
          * @param detector Pointer to the detector for this module instance
          */
         EventLoaderYarr(Configuration& config, std::shared_ptr<Detector> detector);
         ~EventLoaderYarr() {}
 
         /**
          * @brief [Initialise this module] : Look for the data file, open it and initialise histograms
          */
         void initialize() override;
 
         /**
          * @brief [Run the function of this module] : Process a group of trigger events and update the clipboard accordingly
          */
         StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;
 
         /**
          * @brief [Finalise module]
          */
         void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;
 
     private:
         
         // Timing-related member variables:
         int event_number_;                                             // Event number
         bool trigger_tag_timing_;                                      // Toggle for using the tag time
         double previous_time_ = 0;                                     // Previous event time
         int day_offset_ = 0;                                           // Offset for the day rollover
 
         // Input file-related member variables:
         std::string input_directory_;                                  // Input directory with .raw files
         std::shared_ptr<Detector> detector_;                           // Detector pointer
         std::string file_name_;                                        // Name of the raw file
         std::fstream file_handle_;                                     // File stream for binary reading
 
         // Helper struct to hold header information
         struct TriggerHeader {
             uint32_t tag;
             uint16_t l1id;
             uint16_t bcid;
             uint16_t numHits;
             double time;
         };
         
         // Root Plots
         TH2F* hHitMap;                                                  // 2D hitmap for each detector plane
         TH1F* hEventsVsTagTime;                                         // Number of events vs time
         TH1F* hNumHitsVsTagTime;                                        // Number of hits vs time
 

         // ----------------------------------- Helper Functions  -----------------------------------:
 
         // Read the header of the event which contains timing and trigger window information
         TriggerHeader read_header();
         
         // Read the hits (Col, Row, ToT) in an event and creates Pixel objects
         void read_hits(
             PixelVector& pixels,
             double eventTime,
             uint16_t nHits
         );
         
         // Find the raw data files in the input directory
         std::string find_raw_file(
             const std::string& directory, 
             const std::string& detectorName
         );
        
         // Adjust the timestamp to account for day rollover
         double adjust_time(
             const double& current_time
         );
 
     };
 
 } // namespace corryvreckan
 
