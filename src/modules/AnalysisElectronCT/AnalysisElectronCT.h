/**
 * @file
 * @brief Definition of module AnalysisElectronCT
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <iostream>
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
    class AnalysisElectronCT : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        AnalysisElectronCT(Configuration& config, std::shared_ptr<Detector> detector);

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
        std::shared_ptr<Detector> m_detector;
        int m_eventNumber;
        bool charge_weighting_;
        bool fitted_profile_;

        bool ignore_first_frame_;
        bool compute_xray_;
        bool store_data_;

        std::vector<double> widthsXLast1500;
        std::vector<double> widthsYLast1500;

        int x_frames_;
        int chop_frames_;

        // Graphs:
        TH2F* hitMapIntegrated;
        TH2F* hitMapLastFrame;

        TH2F* chargeMapIntegrated;
        TH2F* chargeMapLastFrame;

        TH1F* hitsPerFrame;
        TH1F* chargePerFrame;

        TH1F* projectionX;
        TH1F* projectionY;

        TH1F* projectionChargeX;
        TH1F* projectionChargeY;

        TH1F* widthsX;
        TH1F* widthsY;

        TH1F* centersX;
        TH1F* centersY;

        TProfile* nHitsVsFrame;
        TProfile* chargeVsFrame;

        TH2F* projectionXVsFrame;
        TH2F* projectionYVsFrame;

        TH2F* projectionChargeXVsFrame;
        TH2F* projectionChargeYVsFrame;

        TProfile* centersXVsFrame;
        TProfile* centersYVsFrame;

        TProfile* widthsXVsFrame;
        TProfile* widthsYVsFrame;

        TH1F* widthsXVsFrameLast1500;
        TH1F* widthsYVsFrameLast1500;

        TH1F* framesVsTime;

        TH1F* hitTimeWithinFrame;

        TH1F* pixelCharge;

        TH2F* xray;
    };

} // namespace corryvreckan
