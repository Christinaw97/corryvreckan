/**
 * @file
 * @brief Definition of module Correlations
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRELATIONS_H
#define CORRELATIONS_H 1

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/TimerSignal.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class Correlations : public Module {

    public:
        // Constructors and destructors
        Correlations(Configuration& config, std::shared_ptr<Detector> detector);
        ~Correlations() {}

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        std::shared_ptr<Detector> m_detector;

        // Pixel histograms
        TH2F* hitmap;
        TH2F* hitmap_clusters;
        TH1F* eventTimes;
        TH1F* eventTimesTimerSignal;

        // Correlation plots
        TH1F* correlationX;
        TH1F* correlationXY;
        TH1F* correlationY;
        TH1F* correlationYX;
        TH2F* correlationX2Dlocal;
        TH2F* correlationY2Dlocal;
        TH2F* correlationColCol_px;
        TH2F* correlationColRow_px;
        TH2F* correlationRowCol_px;
        TH2F* correlationRowRow_px;
        TH2F* correlationX2D;
        TH2F* correlationY2D;
        TH2F* correlationYX2D;
        TH2F* correlationXY2D;
        TH2F* correlationXVsTrigger;
        TH2F* correlationYVsTrigger;
        TH2F* correlationYXVsTrigger;
        TH2F* correlationXYVsTrigger;
        TH1F* correlationTime;
        TH2F* correlationTimeOverTime;
        TH2F* correlationTimerSignalTimeOverTime_px;
        TH2F* correlationTimeOverSeedPixelRawValue;
        TH1F* correlationTime_px;
        TH1F* correlationTimerSignalTime_px;
        TH2F* correlationTimeOverTime_px;
        TH2F* correlationTimeOverPixelRawValue_px;
        TH1F* correlationTimeInt;
        TH2F* correlationXVsTime;
        TH2F* correlationYVsTime;
        TH2F* correlationXYVsTime;
        TH2F* correlationYXVsTime;

        // Parameters which can be set by user
        double time_cut_;
        bool do_time_cut_;
        bool corr_vs_time_;
        double time_binning_;

        // Functions for booking of histogram
        // Booking of histograms for normal detectors
        void bookStandardHistograms(int, double, int, std::shared_ptr<Detector>);
        // Booking of histograms for auxiliary detectors
        void bookAuxiliaryHistograms();
    };
} // namespace corryvreckan
#endif // CORRELATIONS_H
