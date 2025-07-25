/**
 * @file
 * @brief Implementation of module Prealignment
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "Prealignment.h"
#include "tools/cuts.h"

using namespace corryvreckan;
using namespace std;

Prealignment::Prealignment(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), m_detector(detector) {

    // Backwards compatibility: also allow timing_cut to be used for time_cut_abs
    config_.setAlias("time_cut_abs", "timing_cut", true);

    config_.setDefault<double>("max_correlation_rms", Units::get<double>(6, "mm"));
    config_.setDefault<double>("damping_factor", 1.0);
    config_.setDefault<PrealignMethod>("method", PrealignMethod::MEAN);
    config_.setDefault<int>("fit_range_rel", 500);
    config_.setDefault<double>("range_abs", Units::get<double>(10, "mm"));
    config_.setDefault<double>("time_range_abs", Units::get<double>(100, "ns"));
    config_.setDefault<double>("time_binning", Units::get<double>(1, "ns"));
    config_.setDefault<int>("nbins_global", 1000);
    config_.setDefault<bool>("align_time", false);

    if(config_.count({"time_cut_rel", "time_cut_abs"}) == 0) {
        config_.setDefault("time_cut_rel", 3.0);
    }

    // timing cut, relative (x * time_resolution) or absolute:
    timeCut = corryvreckan::calculate_cut<double>("time_cut", config_, m_detector);

    max_correlation_rms = config_.get<double>("max_correlation_rms");
    damping_factor = config_.get<double>("damping_factor");
    range_abs = config_.get<double>("range_abs");
    time_range_abs_ = config_.get<double>("time_range_abs");
    nbins_global = config_.get<int>("nbins_global");
    method = config_.get<PrealignMethod>("method");
    fit_range_rel = config_.get<int>("fit_range_rel");
    fixed_planes_ = config_.getArray<std::string>("fixed_planes", {});
    align_time_ = config_.get<bool>("align_time");
    time_binning_ = config_.get<double>("time_binning");

    LOG(DEBUG) << "Setting max_correlation_rms to : " << max_correlation_rms;
    LOG(DEBUG) << "Setting damping_factor to : " << damping_factor;
}

void Prealignment::initialize() {

    LOG(INFO) << "Spatial histograms have a range of " << Units::display(range_abs, {"um", "mm"}) << " and " << nbins_global
              << " bins";
    LOG(INFO) << "Time histogram has a range of " << Units::display(time_range_abs_, {"ps", "ns", "us"}) << " and binning "
              << Units::display(time_binning_, {"ps", "ns", "us"});

    // get the reference detector:
    std::shared_ptr<Detector> reference = get_reference();

    // Correlation plots
    std::string title = m_detector->getName() + ": correlation X;x_{ref}-x [mm];events";
    correlationX = new TH1F("correlationX", title.c_str(), nbins_global, -1.0 * range_abs, 1.0 * range_abs);
    title = m_detector->getName() + ": correlation Y;y_{ref}-y [mm];events";
    correlationY = new TH1F("correlationY", title.c_str(), nbins_global, -1.0 * range_abs, 1.0 * range_abs);
    title = m_detector->getName() + ": correlation XY;x_{ref}-x [mm];y_{ref}-y [mm];events";
    correlationXY = new TH2F("correlationXY",
                             title.c_str(),
                             nbins_global,
                             -1.0 * range_abs,
                             1.0 * range_abs,
                             nbins_global,
                             -1.0 * range_abs,
                             1.0 * range_abs);
    if(align_time_) {
        title = m_detector->getName() + ": correlation time;t_{ref}-t [ns];events";
        const auto time_bins = static_cast<int>(time_range_abs_ / time_binning_);
        correlationTime_ =
            new TH1F("correlationTime", title.c_str(), time_bins, -1.0 * time_range_abs_, 1.0 * time_range_abs_);
    }
    // 2D correlation plots (pixel-by-pixel, local coordinates):
    title = m_detector->getName() + ": 2D correlation X (local);x [px];x_{ref} [px];events";
    correlationX2Dlocal = new TH2F("correlationX_2Dlocal",
                                   title.c_str(),
                                   m_detector->nPixels().X(),
                                   -0.5,
                                   m_detector->nPixels().X() - 0.5,
                                   reference->nPixels().X(),
                                   -0.5,
                                   reference->nPixels().X() - 0.5);
    title = m_detector->getName() + ": 2D correlation Y (local);y [px];y_{ref} [px];events";
    correlationY2Dlocal = new TH2F("correlationY_2Dlocal",
                                   title.c_str(),
                                   m_detector->nPixels().Y(),
                                   -0.5,
                                   m_detector->nPixels().Y() - 0.5,
                                   reference->nPixels().Y(),
                                   -0.5,
                                   reference->nPixels().Y() - 0.5);
    title = m_detector->getName() + ": 2D correlation X (global);x [mm];x_{ref} [mm];events";
    correlationX2D = new TH2F("correlationX_2D",
                              title.c_str(),
                              nbins_global / 10,
                              -1.0 * range_abs,
                              1.0 * range_abs,
                              nbins_global / 10,
                              -1.0 * range_abs,
                              1.0 * range_abs);
    title = m_detector->getName() + ": 2D correlation Y (global);y [mm];y_{ref} [mm];events";
    correlationY2D = new TH2F("correlationY_2D",
                              title.c_str(),
                              nbins_global / 10,
                              -1.0 * range_abs,
                              1.0 * range_abs,
                              nbins_global / 10,
                              -1.0 * range_abs,
                              1.0 * range_abs);
}

StatusCode Prealignment::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Get the clusters
    auto clusters = clipboard->getData<Cluster>(m_detector->getName());
    if(clusters.empty()) {
        LOG(DEBUG) << "Detector " << m_detector->getName() << " does not have any clusters on the clipboard";
        return StatusCode::NoData;
    }

    // Get clusters from reference detector
    auto reference = get_reference();
    auto referenceClusters = clipboard->getData<Cluster>(reference->getName());
    if(referenceClusters.empty()) {
        LOG(DEBUG) << "Reference detector " << reference->getName() << " does not have any clusters on the clipboard";
        return StatusCode::NoData;
    }

    // Loop over all clusters and fill histograms
    for(auto& cluster : clusters) {
        // Loop over reference plane pixels to make correlation plots
        for(auto& refCluster : referenceClusters) {
            double timeDifference = refCluster->timestamp() - cluster->timestamp();

            // Correlation plots
            if(abs(timeDifference) < timeCut) {
                correlationX->Fill(refCluster->global().x() - cluster->global().x());
                correlationX2D->Fill(cluster->global().x(), refCluster->global().x());
                correlationX2Dlocal->Fill(cluster->column(), refCluster->column());
                correlationY->Fill(refCluster->global().y() - cluster->global().y());
                correlationY2D->Fill(cluster->global().y(), refCluster->global().y());
                correlationY2Dlocal->Fill(cluster->row(), refCluster->row());
                correlationXY->Fill(refCluster->global().x() - cluster->global().x(),
                                    refCluster->global().y() - cluster->global().y());
                if(align_time_) {
                    correlationTime_->Fill(refCluster->timestamp() - cluster->timestamp());
                }
            }
        }
    }

    return StatusCode::Success;
}

void Prealignment::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    double rmsX = correlationX->GetRMS();
    double rmsY = correlationY->GetRMS();
    if(rmsX > max_correlation_rms or rmsY > max_correlation_rms) {
        LOG(ERROR) << "Detector " << m_detector->getName() << ": RMS is too wide for prealignment shifts";
        LOG(ERROR) << "Detector " << m_detector->getName() << ": RMS X = " << Units::display(rmsX, {"mm", "um"})
                   << " , RMS Y = " << Units::display(rmsY, {"mm", "um"});
    }
    bool is_fixed = std::find(fixed_planes_.begin(), fixed_planes_.end(), m_detector->getName()) != fixed_planes_.end();

    // Move all but the reference and user-defined plane:
    if(!m_detector->isReference() && !is_fixed) {
        LOG(INFO) << "Running detector " << m_detector->getName();

        double shift_X = 0.;
        double shift_Y = 0.;
        double shift_T = 0.;

        LOG(INFO) << "Using prealignment method: " << corryvreckan::to_string(method);
        if(method == PrealignMethod::GAUSS_FIT) {
            int binMaxX = correlationX->GetMaximumBin();
            double fit_low_x =
                correlationX->GetXaxis()->GetBinCenter(binMaxX) - m_detector->getSpatialResolution().x() * fit_range_rel;
            double fit_high_x =
                correlationX->GetXaxis()->GetBinCenter(binMaxX) + m_detector->getSpatialResolution().x() * fit_range_rel;

            int binMaxY = correlationY->GetMaximumBin();
            double fit_low_y =
                correlationY->GetXaxis()->GetBinCenter(binMaxY) - m_detector->getSpatialResolution().y() * fit_range_rel;
            double fit_high_y =
                correlationY->GetXaxis()->GetBinCenter(binMaxY) + m_detector->getSpatialResolution().y() * fit_range_rel;

            LOG(DEBUG) << "Fit range in x direction from: " << Units::display(fit_low_x, {"mm", "um"}) << " to "
                       << Units::display(fit_high_x, {"mm", "um"});
            LOG(DEBUG) << "Fit range in y direction from: " << Units::display(fit_low_y, {"mm", "um"}) << " to "
                       << Units::display(fit_high_y, {"mm", "um"});

            correlationX->Fit("gaus", "Q", "", fit_low_x, fit_high_x);
            correlationY->Fit("gaus", "Q", "", fit_low_y, fit_high_y);

            shift_X = correlationX->GetFunction("gaus")->GetParameter(1);
            shift_Y = correlationY->GetFunction("gaus")->GetParameter(1);

            if(align_time_) {
                int binMaxTime = correlationTime_->GetMaximumBin();
                double fit_low_t =
                    correlationTime_->GetXaxis()->GetBinCenter(binMaxTime) - m_detector->getTimeResolution() * fit_range_rel;
                double fit_high_t =
                    correlationTime_->GetXaxis()->GetBinCenter(binMaxTime) + m_detector->getTimeResolution() * fit_range_rel;
                correlationTime_->Fit("gaus", "Q", "", fit_low_t, fit_high_t);
                shift_T = correlationTime_->GetFunction("gaus")->GetParameter(1);
            }
        } else if(method == PrealignMethod::MEAN) {
            shift_X = correlationX->GetMean();
            shift_Y = correlationY->GetMean();
            if(align_time_) {
                shift_T = correlationTime_->GetMean();
            }
        } else if(method == PrealignMethod::MAXIMUM) {
            int binMaxX = correlationX->GetMaximumBin();
            shift_X = correlationX->GetXaxis()->GetBinCenter(binMaxX);
            int binMaxY = correlationY->GetMaximumBin();
            shift_Y = correlationY->GetXaxis()->GetBinCenter(binMaxY);
            if(align_time_) {
                int binMaxTime = correlationTime_->GetMaximumBin();
                shift_T = correlationTime_->GetXaxis()->GetBinCenter(binMaxTime);
            }
        } else if(method == PrealignMethod::MAXIMUM2D) {
            int binMaxX1 = correlationX->GetMaximumBin();
            TH1D* ProjY = correlationXY->ProjectionY("_py", binMaxX1 - 1, binMaxX1 + 1);
            int binMaxY1 = ProjY->GetMaximumBin();
            auto max1 = correlationXY->GetBinContent(binMaxX1, binMaxY1);

            int binMaxY2 = correlationY->GetMaximumBin();
            TH1D* ProjX = correlationXY->ProjectionX("_px", binMaxY2 - 1, binMaxY2 + 1);
            int binMaxX2 = ProjX->GetMaximumBin();
            auto max2 = correlationXY->GetBinContent(binMaxX2, binMaxY2);

            if(max1 > max2) {
                shift_X = correlationX->GetXaxis()->GetBinCenter(binMaxX1);
                shift_Y = ProjY->GetXaxis()->GetBinCenter(binMaxY1);
            } else {
                shift_X = ProjX->GetXaxis()->GetBinCenter(binMaxX2);
                shift_Y = correlationY->GetXaxis()->GetBinCenter(binMaxY2);
            }
        }

        LOG(DEBUG) << "Shift (without damping factor)" << m_detector->getName()
                   << ": x = " << Units::display(shift_X, {"mm", "um"}) << " , y = " << Units::display(shift_Y, {"mm", "um"})
                   << " , t = " << Units::display(shift_T, {"ns"});
        LOG(INFO) << "Move in x by = " << Units::display(shift_X * damping_factor, {"mm", "um"})
                  << " , and in y by = " << Units::display(shift_Y * damping_factor, {"mm", "um"})
                  << " , and in t by = " << Units::display(shift_T * damping_factor, {"ns"});
        LOG(INFO) << "Detector position after shift in x = "
                  << Units::display(m_detector->displacement().X() + damping_factor * shift_X, {"mm", "um"})
                  << " , and in y = "
                  << Units::display(m_detector->displacement().Y() + damping_factor * shift_Y, {"mm", "um"})
                  << " , and in t = " << Units::display(m_detector->timeOffset() + damping_factor * shift_T, {"ns"});
        m_detector->update(XYZPoint(m_detector->displacement().X() + damping_factor * shift_X,
                                    m_detector->displacement().Y() + damping_factor * shift_Y,
                                    m_detector->displacement().Z()),
                           m_detector->rotation());
        m_detector->setTimeOffset(m_detector->timeOffset() + damping_factor * shift_T);
    }
}
