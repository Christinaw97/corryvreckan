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
#include <TStyle.h>
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
        int m_eventNumber;
        uint64_t m_absolute_event_time;

        // Input file-related member variables:
        std::string m_inputDirectory;               // Input directory with .raw files
        std::shared_ptr<Detector> m_detector;       // Detector pointer
        std::string m_filename;                      // Name of the raw file
        std::fstream fileHandle;                    // File stream for binary reading

        uint64_t m_previous_time = 0;               // Previous event time
        uint64_t m_day_offset = 0;                  // Offset for the day rollover

        // Helper struct to hold header information
        struct TriggerHeader {
            uint32_t tag;
            uint16_t l1id;
            uint16_t bcid;
            uint16_t numHits;
            uint8_t basetag;
            uint8_t exttag;
            uint32_t time;
        };
        
        // Root Plots
        TH2F* hHitMap;                              // 2D hitmap for each detector plane
        TH1F* numHitsVsTime;                        // Number of hits vs time
        TH1F* eventsVsAbsoluteTime;                 // Number of events vs absolute time
        TH1F* numHitsVsAbsTime;                     // Number of hits vs absolute time

        // Helper Functions:

        // Read the header of the event
        TriggerHeader readHeader();
        
        // Read the hits of the event
        void readHits(
            PixelVector& pixels,
            uint32_t eventTime,
            uint16_t nHits
        );
        
        // Find the raw file in the input directory
        std::string findRawFile(
            const std::string& directory, 
            const std::string& detectorName
        );

        // Compute the absolute event time
        uint64_t adjustEventTime(
            const std::vector<uint32_t>& triggerTimes
        );

    };

} // namespace corryvreckan
