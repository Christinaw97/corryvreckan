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
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class EventLoaderYarr : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        EventLoaderYarr(Configuration& config, std::shared_ptr<Detector> detector);
        ~EventLoaderYarr() {}

        /**
         * @brief [Initialise this module]
         */
        void initialize() override;

        /**
         * @brief [Run the function of this module]
         */
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

        /**
         * @brief [Finalise module]
         */
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

    private:
        int m_eventNumber;

        // Event buffer
        uint32_t this_tag, this_time; 
        uint16_t this_l1id, this_bcid, this_t_hits;
        uint8_t this_basetag, this_exttag;

        // Member variables
        std::string m_inputDirectory;           // Input directory
        std::shared_ptr<Detector> m_detector;   // Detector
        std::string filename_;                  // Current file name
        std::streampos filePos;                 // File position
        std::fstream fileHandle;                // File handle
        int m_triggerMultiplicity;              // Trigger multiplicity
        uint64_t m_bufferDepth;                 // Buffer depth
        bool header_read;                       // Header read flag
        
        // Root Plots
        TH2F* hHitMap;                          // Make a hitmap for each event

        // Helper Functions
        void readHeader();                      // Read the header of the event
        void readHits(PixelVector& pixels);     // Read the hits of the event

    };

} // namespace corryvreckan
