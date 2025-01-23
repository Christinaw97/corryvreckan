/**
 * @file
 * @brief Implementation of [AnalysisItkStripEfficiency] module

 */

#include "AnalysisItkStripEfficiency.h"

#include "core/detector/PolarDetector.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

using namespace corryvreckan;

AnalysisItkStripEfficiency::AnalysisItkStripEfficiency(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector) {
    m_detector = detector;

    config_.setDefault<double>("time_cut_frameedge", Units::get<double>(20, "ns"));
    config_.setDefault<double>("chi2ndof_cut", 3.);
    config_.setDefault<int>("profile_bins", 300);
    config_.setDefault<double>("profile_xrange", 1.5 * m_detector->getSize().X());
    config_.setDefault<double>("profile_yrange", 1.5 * m_detector->getSize().Y());
    config_.setDefault<double>("perimeter_exclude", 1.);
    config_.setDefault<double>("inpixel_bin_size", Units::get<double>(1.0, "um"));
    config_.setDefault<XYVector>("inpixel_cut_edge", {Units::get(5.0, "um"), Units::get(5.0, "um")});
    config_.setDefault<double>("masked_pixel_distance_cut", 1.);
    config_.setDefault<std::string>("file_ttc", "");
    config_.setDefault<std::string>("ttc_tag", "PTDC_DUT.BIT");
    config_.setDefault<std::string>("eudaq_loglevel", "ERROR");
    // config_.setDefault<std::vector<double>>("delay_cut", {Units::get<double>(14, "ns"), Units::get<double>(19, "ns")});
    std::vector<int> defaultLimits = {0, 64};
    config_.setDefaultArray<int>("delay_cuts", defaultLimits);

    readerFile_ = config_.getPath("file_ttc");
    ttc_tag_ = config_.get<std::string>("ttc_tag");
    m_timeCutFrameEdge = config_.get<double>("time_cut_frameedge");
    m_chi2ndofCut = config_.get<double>("chi2ndof_cut");
    m_perimeter_exclude = config_.get<double>("perimeter_exclude");
    m_inpixelBinSize = config_.get<double>("inpixel_bin_size");
    require_associated_cluster_on_ = config_.getArray<std::string>("require_associated_cluster_on", {});
    m_inpixelEdgeCut = config_.get<XYVector>("inpixel_cut_edge");
    m_maskedPixelDistanceCut = config_.get<int>("masked_pixel_distance_cut");
    // m_delay_cut = config_.get<std::vector<double>>("delay_cut");
    m_delay_cuts = config_.getArray<int>("delay_cuts");
    m_profile_bins = config_.get<int>("profile_bins");
    m_profile_xrange = config_.get<double>("profile_xrange");
    m_profile_yrange = config_.get<double>("profile_yrange");

    LOG(INFO) << "time_cut_frameedge = " << m_timeCutFrameEdge;

    LOG(DEBUG) << "m_detector->getSize().X() = " << m_detector->getSize().X();
    LOG(DEBUG) << "m_detector->getSize().Y() = " << m_detector->getSize().Y();

    if(m_delay_cuts.size() == 1) {
        m_delay_cuts.push_back(0);
    } else {
        if(m_delay_cuts.size() > 2) {
            LOG(INFO) << "More than 2 values read for delay limits, using 2 only";
        }
    }

    // Set EUDAQ log level to desired value:
    EUDAQ_LOG_LEVEL(config_.get<std::string>("eudaq_loglevel"));
    LOG(INFO) << "Setting EUDAQ2 log level to \"" << config_.get<std::string>("eudaq_loglevel") << "\"";

    // Prepare EUDAQ2 config object
    eudaq::Configuration cfg;

    // Forward all settings to EUDAQ
    // WARNING: the EUDAQ Configuration class is not very flexible and e.g. booleans have to be passed as 1 and 0.
    auto configs = config_.getAll();
    for(const auto& key : configs) {
        LOG(DEBUG) << "Forwarding key \"" << key.first << " = " << key.second << "\" to EUDAQ converter";
        cfg.Set(key.first, key.second);
    }

    // Converting the newly built configuration to a shared pointer of a const configuration object
    // Unfortunately, EUDAQ does not provide appropriate member functions for their configuration class to avoid this dance
    const eudaq::Configuration eu_cfg = cfg;
    eudaq_config_ = std::make_shared<const eudaq::Configuration>(eu_cfg);
}

void AnalysisItkStripEfficiency::initialize() {

    LOG(INFO) << "Opening TTC stream file: " << readerFile_;
    try {
        readerTTC_ = eudaq::Factory<eudaq::FileReader>::MakeUnique(eudaq::str2hash("native"), readerFile_);
    } catch(...) {
        LOG(ERROR) << "eudaq::FileReader could not read the input file ' " << readerFile_
                   << " '. Please verify that the path and file name are correct.";
        throw InvalidValueError(config_, "file_path", "Parsing error!");
    }

    LOG(INFO) << "Setting track delay window to";
    LOG(INFO) << "Min delay cut is: " << m_delay_cuts[1];
    LOG(INFO) << "Max delay cut is: " << m_delay_cuts[0];

    eTotalEfficiency = new TEfficiency("eTotalEfficiency", "totalEfficiency;;#epsilon", 1, 0, 1);
    eTotalEfficiency->SetDirectory(this->getROOTDirectory());
    // LOG(DEBUG) << "eTotalEfficiency pointer print" << eTotalEfficiency;

    eTimingEfficiency = new TEfficiency("eTimingEfficiency", "TimingEfficiency;Delay;#epsilon", 52, -0.5, 51.5);
    eTimingEfficiency->SetDirectory(this->getROOTDirectory());
    // LOG(DEBUG) << "eTimingEfficiency pointer print" << eTimingEfficiency;

    eTotalEfficiency_inPixelROI = new TEfficiency(
        "eTotalEfficiency_inPixelROI", "eTotalEfficiency_inPixelROI;;#epsilon (within in-pixel ROI)", 1, 0, 1);

    // LOG(DEBUG) << "eTotalEfficiency_inPixelROI pointer print" << eTotalEfficiency_inPixelROI;

    hPixelEfficiency = new TH1D(
        "hPixelEfficiency", "hPixelEfficiency; single pixel efficiency; # entries", 201, 0, 1.005); // get 0.5%-wide bins
    hPixelEfficiency->SetDirectory(this->getROOTDirectory());

    hPixelEfficiencyMatrix = new TH1D("hPixelEfficiencyMatrix",
                                      "hPixelEfficiencyMatrix; single pixel efficiency; # entries",
                                      201,
                                      0,
                                      1.005); // get 0.5%-wide bins
    hPixelEfficiencyMatrix->SetDirectory(this->getROOTDirectory());

    pitch_x = static_cast<double>(Units::convert(m_detector->getPitch().X(), "mrad"));
    pitch_y = static_cast<double>(Units::convert(m_detector->getPitch().Y(), "um"));

    auto nbins_x = static_cast<int>(std::ceil(m_detector->getPitch().X() / m_inpixelBinSize));
    auto nbins_y = static_cast<int>(std::ceil(m_detector->getPitch().Y() / m_inpixelBinSize));
    if(nbins_x > 1e4 || nbins_y > 1e4) {
        throw InvalidValueError(config_, "inpixel_bin_size", "Too many bins for in-pixel histograms.");
    }

    auto title = m_detector->getName() + " Pixel efficiency map;in-pixel x_{track} [#mum];in-pixel y_{track} #mum;#epsilon";
    hPixelEfficiencyMap_trackPos_TProfile = new TProfile2D("pixelEfficiencyMap_trackPos_TProfile",
                                                           title.c_str(),
                                                           nbins_x,
                                                           -pitch_x / 2.,
                                                           pitch_x / 2.,
                                                           nbins_y,
                                                           -pitch_y / 2.,
                                                           pitch_y / 2.,
                                                           0,
                                                           1);

    hPixelEfficiencyMap_trackPos = new TEfficiency("pixelEfficiencyMap_trackPos",
                                                   title.c_str(),
                                                   nbins_x,
                                                   -pitch_x / 2.,
                                                   pitch_x / 2.,
                                                   nbins_y,
                                                   -pitch_y / 2.,
                                                   pitch_y / 2.);
    hPixelEfficiencyMap_trackPos->SetDirectory(this->getROOTDirectory());

    LOG(DEBUG) << "hPixelEfficiencyMap_trackPos point print" << hPixelEfficiencyMap_trackPos;

    title = m_detector->getName() +
            " Pixel efficiency map (in-pixel ROI);in-pixel x_{track} [#mum];in-pixel y_{track} #mum;#epsilon";
    hPixelEfficiencyMap_inPixelROI_trackPos_TProfile = new TProfile2D("pixelEfficiencyMap_inPixelROI_trackPos_TProfile",
                                                                      title.c_str(),
                                                                      nbins_x,
                                                                      -pitch_x / 2.,
                                                                      pitch_x / 2.,
                                                                      nbins_y,
                                                                      -pitch_y / 2.,
                                                                      pitch_y / 2.,
                                                                      0,
                                                                      1);

    title = m_detector->getName() + " Chip efficiency map;x [px];y [px];#epsilon";
    hChipEfficiencyMap_trackPos_TProfile = new TProfile2D("chipEfficiencyMap_trackPos_TProfile",
                                                          title.c_str(),
                                                          m_detector->nPixels().X(),
                                                          -0.5,
                                                          m_detector->nPixels().X() - 0.5,
                                                          m_detector->nPixels().Y(),
                                                          -0.5,
                                                          m_detector->nPixels().Y() - 0.5,
                                                          0,
                                                          1);

    hChipEfficiencyMap_trackPos = new TEfficiency("chipEfficiencyMap_trackPos",
                                                  title.c_str(),
                                                  m_detector->nPixels().X(),
                                                  -0.5,
                                                  m_detector->nPixels().X() - 0.5,
                                                  m_detector->nPixels().Y(),
                                                  -0.5,
                                                  m_detector->nPixels().Y() - 0.5);
    hChipEfficiencyMap_trackPos->SetDirectory(this->getROOTDirectory());
    LOG(DEBUG) << "hChipEfficiencyMap_trackPos point print" << hChipEfficiencyMap_trackPos;

    title = m_detector->getName() + " Pixel efficiency matrix;x [px];y [px];#epsilon";
    hPixelEfficiencyMatrix_TProfile = new TProfile2D("hPixelEfficiencyMatrixTProfile",
                                                     title.c_str(),
                                                     m_detector->nPixels().X(),
                                                     -0.5,
                                                     m_detector->nPixels().X() - 0.5,
                                                     m_detector->nPixels().Y(),
                                                     -0.5,
                                                     m_detector->nPixels().Y() - 0.5,
                                                     0,
                                                     1);

    title = m_detector->getName() + " Global efficiency map;x [mm];y [mm];#epsilon";
    hGlobalEfficiencyMap_trackPos_TProfile = new TProfile2D("globalEfficiencyMap_trackPos_TProfile",
                                                            title.c_str(),
                                                            m_profile_bins,
                                                            -m_profile_xrange,
                                                            m_profile_xrange,
                                                            m_profile_bins,
                                                            -m_profile_yrange,
                                                            m_profile_yrange,
                                                            0,
                                                            1);
    hGlobalEfficiencyMap_trackPos_TProfile->SetDirectory(this->getROOTDirectory());

    LOG(DEBUG) << " m_profile_bins: " << m_profile_bins;
    LOG(DEBUG) << " m_profile_xrange: " << m_profile_xrange;
    LOG(DEBUG) << "-m_profile_xrange: " << -m_profile_xrange;
    LOG(DEBUG) << "-m_profile_yrange: " << -m_profile_yrange;
    LOG(DEBUG) << " m_profile_yrange: " << m_profile_yrange;
    hGlobalEfficiencyMap_trackPos = new TEfficiency("globalEfficiencyMap_trackPos",
                                                    title.c_str(),
                                                    m_profile_bins,
                                                    -m_profile_xrange,
                                                    m_profile_xrange,
                                                    m_profile_bins,
                                                    -m_profile_yrange,
                                                    m_profile_yrange);
    hGlobalEfficiencyMap_trackPos->SetDirectory(this->getROOTDirectory());

    title = m_detector->getName() + " Chip efficiency map;x [px];y [px];#epsilon";
    hChipEfficiencyMap_clustPos_TProfile = new TProfile2D("chipEfficiencyMap_clustPos_TProfile",
                                                          title.c_str(),
                                                          m_detector->nPixels().X(),
                                                          -0.5,
                                                          m_detector->nPixels().X() - 0.5,
                                                          m_detector->nPixels().Y(),
                                                          -0.5,
                                                          m_detector->nPixels().Y() - 0.5,
                                                          0,
                                                          1);
    hChipEfficiencyMap_clustPos = new TEfficiency("chipEfficiencyMap_clustPos",
                                                  title.c_str(),
                                                  m_detector->nPixels().X(),
                                                  -0.5,
                                                  m_detector->nPixels().X() - 0.5,
                                                  m_detector->nPixels().Y(),
                                                  -0.5,
                                                  m_detector->nPixels().Y() - 0.5);
    hChipEfficiencyMap_clustPos->SetDirectory(this->getROOTDirectory());

    title = m_detector->getName() + " Global efficiency map;x [mm];y [mm];#epsilon";
    hGlobalEfficiencyMap_clustPos_TProfile = new TProfile2D("globalEfficiencyMap_clustPos_TProfile",
                                                            title.c_str(),
                                                            m_profile_bins,
                                                            -m_profile_xrange,
                                                            m_profile_xrange,
                                                            m_profile_bins,
                                                            -m_profile_yrange,
                                                            m_profile_yrange,
                                                            0,
                                                            1);

    hGlobalEfficiencyMap_clustPos = new TEfficiency("globalEfficiencyMap_clustPos",
                                                    title.c_str(),
                                                    m_profile_bins,
                                                    -m_profile_xrange,
                                                    m_profile_xrange,
                                                    m_profile_bins,
                                                    -m_profile_yrange,
                                                    m_profile_yrange);
    hGlobalEfficiencyMap_clustPos->SetDirectory(this->getROOTDirectory());

    hDistanceCluster = new TH1D("distanceTrackHit",
                                "distance between track and hit; | #vec{track} - #vec{dut} | [mm]",
                                static_cast<int>(std::sqrt(m_detector->getPitch().x() * m_detector->getPitch().y())),
                                0,
                                std::sqrt(m_detector->getPitch().x() * m_detector->getPitch().y()));
    hDistanceCluster_track = new TH2D("distanceTrackHit2D",
                                      "distance between track and hit; track_x - dut_x [mm]; track_y - dut_y [mm] ",
                                      150,
                                      -1.5 * m_detector->getPitch().x(),
                                      1.5 * m_detector->getPitch().x(),
                                      150,
                                      -1.5 * m_detector->getPitch().y(),
                                      1.5 * m_detector->getPitch().y());

    efficiencyColumns = new TEfficiency("efficiencyColumns",
                                        "Efficiency vs. column number; column; #epsilon",
                                        m_detector->nPixels().X(),
                                        -0.5,
                                        m_detector->nPixels().X() - 0.5);
    efficiencyColumns->SetDirectory(this->getROOTDirectory());

    efficiencyRows = new TEfficiency("efficiencyRows",
                                     "Efficiency vs. row number; row; #epsilon",
                                     m_detector->nPixels().Y(),
                                     -0.5,
                                     m_detector->nPixels().Y() - 0.5);
    efficiencyRows->SetDirectory(this->getROOTDirectory());

    efficiencyVsTime = new TEfficiency("efficiencyVsTime", "Efficiency vs. time; time [s]; #epsilon", 3000, 0, 3000);
    efficiencyVsTime->SetDirectory(this->getROOTDirectory());

    hTrackTimeToPrevHit_matched =
        new TH1D("trackTimeToPrevHit_matched", "trackTimeToPrevHit_matched;time to prev hit [us];# events", 1e6, 0, 1e6);
    hTrackTimeToPrevHit_notmatched = new TH1D(
        "trackTimeToPrevHit_notmatched", "trackTimeToPrevHit_notmatched;time to prev hit [us];# events", 1e6, 0, 1e6);

    title = m_detector->getName() + "time difference to previous track (if this has assoc cluster)";
    hTimeDiffPrevTrack_assocCluster = new TH1D("timeDiffPrevTrack_assocCluster", title.c_str(), 11000, -1000, 10000);
    hTimeDiffPrevTrack_assocCluster->GetXaxis()->SetTitle("time diff [#mus]");
    hTimeDiffPrevTrack_assocCluster->GetYaxis()->SetTitle("events");
    title = m_detector->getName() + "time difference to previous track (if this has no assoc cluster)";
    hTimeDiffPrevTrack_noAssocCluster = new TH1D("timeDiffPrevTrack_noAssocCluster", title.c_str(), 11000, -1000, 10000);
    hTimeDiffPrevTrack_noAssocCluster->GetXaxis()->SetTitle("time diff [#mus]");
    hTimeDiffPrevTrack_noAssocCluster->GetYaxis()->SetTitle("events");

    hRowDiffPrevTrack_assocCluster =
        new TH1D("rowDiffPrevTrack_assocCluster",
                 "rowDiffPrevTrack_assocCluster; row difference (matched track to prev track) [px];# events",
                 2 * m_detector->nPixels().Y(),
                 -m_detector->nPixels().Y() - 0.5,
                 m_detector->nPixels().Y() - 0.5);

    hColDiffPrevTrack_assocCluster =
        new TH1D("colDiffPrevTrack_assocCluster",
                 "colDiffPrevTrack_assocCluster;column difference (matched track to prev track) [px];# events",
                 2 * m_detector->nPixels().X(),
                 -m_detector->nPixels().X() - 0.5,
                 m_detector->nPixels().X() - 0.5);

    hRowDiffPrevTrack_noAssocCluster =
        new TH1D("rowDiffPrevTrack_noAssocCluster",
                 "rowDiffPrevTrack_noAssocCluster;row difference (non-matched track - prev track) [px];# events",
                 2 * m_detector->nPixels().Y(),
                 -m_detector->nPixels().Y() - 0.5,
                 m_detector->nPixels().Y() - 0.5);

    hColDiffPrevTrack_noAssocCluster =
        new TH1D("colDiffPrevTrack_noAassocCluster",
                 "colDiffPrevTrack_noAssocCluster;column difference (non-matched track - prev track) [px];# events",
                 2 * m_detector->nPixels().X(),
                 -m_detector->nPixels().X() - 0.5,
                 m_detector->nPixels().X() - 0.5);

    hPosDiffPrevTrack_assocCluster = new TH2D("posDiffPrevTrack_assocCluster",
                                              "posDiffPrevTrack_assocCluster;column difference (matched track - prev track) "
                                              "[px];row difference (matched track - prev track) [px];# events",
                                              2 * m_detector->nPixels().X(),
                                              -m_detector->nPixels().X() - 0.5,
                                              m_detector->nPixels().X() - 0.5,
                                              2 * m_detector->nPixels().Y(),
                                              -m_detector->nPixels().Y() - 0.5,
                                              m_detector->nPixels().Y() - 0.5);

    hPosDiffPrevTrack_noAssocCluster =
        new TH2D("posDiffPrevTrack_noAssocCluster",
                 "posDiffPrevTrack_noAssocCluster;column difference (non-matched track - prev track) [px];row difference "
                 "(non-matched track - prev track) [px];# events",
                 2 * m_detector->nPixels().X(),
                 -m_detector->nPixels().X() - 0.5,
                 m_detector->nPixels().X() - 0.5,
                 2 * m_detector->nPixels().Y(),
                 -m_detector->nPixels().Y() - 0.5,
                 m_detector->nPixels().Y() - 0.5);

    hPos_TrackLocal_Ass = new TH2D("pos_TrackLocal_Ass",
                                   "Local Track position with associated cluster;Track X;Track Y;# events",
                                   200,
                                   -10,
                                   10,
                                   200,
                                   20,
                                   50);
    hPos_TrackLocal_No_Ass = new TH2D("pos_TrackLocal_No_Ass",
                                      "Local Track position without associated cluster;Track X;Track Y;# events",
                                      200,
                                      -10,
                                      10,
                                      200,
                                      20,
                                      50);

    // hresi_x_event_Start = new TH2D("resi_x_event_Start",
    //                                   "X residual vs event start timestamp;Event Time [s];X residual [mrad];# events",
    //                                   1206,
    //                                   -6,
    //                                   600,
    //                                   600,
    //                                   -3.,
    //                                   3.);

    // hresi_x_event_end = new TH2D("resi_x_event_end",
    //                                   "X residual vs event end timestamp;Event Time [s];X residual [mrad];# events",
    //                                   1206,
    //                                   -6,
    //                                   600,
    //                                   600,
    //                                   -3,
    //                                   3);

    // hresi_x_track_time = new TH2D("resi_x_track_time",
    //                                   "X residual vs track timestamp (unit bugged,);Event Time [ns];X residual [mrad];#
    //                                   events", 2100, -100, 600, 600, -3, 3);

    // if (m_detector->getCoordinateType() == "polar") {
    // auto m_detector_polar = std::static_pointer_cast<ITkStripR0>(m_detector);
    // ai to resolve at telescope resolution 3-5 micron. ITk EC strip pitch ~ 70microns, 72/3 = 24
    // for 2 strips plot 2*24 = 48
    title = "InPixelEfficiency;In-2-pixel-position;#epsilon";
    eInPixelEfficiency = new TEfficiency("eInPixelEfficiency", title.c_str(), 48, 0.0, 2.0);
    eInPixelEfficiency->SetDirectory(this->getROOTDirectory());
    // }
    // else{
    //     title = "InPixelEfficiency;In-2-pixel-position;#epsilon";
    //     eInPixelEfficiency =
    //         new TEfficiency("eInPixelEfficiency",title.c_str(), 100, 0.0, 2.0);
    //     eInPixelEfficiency->SetDirectory(this->getROOTDirectory());
    // }

    title = m_detector->getName() + " even-odd efficiency;x [px];y [px];#epsilon";
    hStripEfficiencyOddEven_TProfile =
        new TProfile2D("hStripEfficiencyOddEven", title.c_str(), 2, -0.5, 1.5, 1000, 0, m_detector->getPitch().Y(), 0, 1);
    // initialize matrix with hit timestamps to all 0:
    auto nRows = static_cast<size_t>(m_detector->nPixels().Y());
    auto nCols = static_cast<size_t>(m_detector->nPixels().X());
    std::vector<double> v_row(nRows, 0.); // create vector will zeros of length <nRows>
    prev_hit_ts.assign(nCols, v_row);     // use vector v_row to construct matrix
}

int AnalysisItkStripEfficiency::get_highest_bit_set(int bitset) {
    int index = 0;
    while((1 << index++) & bitset)
        ;
    return index - 1;
}

Event::Position AnalysisItkStripEfficiency::is_within_event(const std::shared_ptr<Clipboard>& clipboard,
                                                            std::shared_ptr<eudaq::StandardEvent> evt) const {
    // Potentially shift the trigger IDs if requested
    auto triggerN = static_cast<uint32_t>(static_cast<int>(evt->GetTriggerN()));

    auto trigger_position = clipboard->getEvent()->getTriggerPosition(triggerN);
    if(trigger_position == Event::Position::BEFORE) {
        LOG(DEBUG) << "Trigger ID " << evt->GetTriggerN() << " before triggers registered in Corryvreckan event";
        LOG(DEBUG) << "(Shifted) Trigger ID " << triggerN << " before triggers registered in Corryvreckan event";
    } else if(trigger_position == Event::Position::AFTER) {
        LOG(DEBUG) << "Trigger ID " << evt->GetTriggerN() << " after triggers registered in Corryvreckan event";
        LOG(DEBUG) << "(Shifted) Trigger ID " << triggerN << " after triggers registered in Corryvreckan event";
    } else if(trigger_position == Event::Position::UNKNOWN) {
        LOG(DEBUG) << "Trigger ID " << evt->GetTriggerN() << " within Corryvreckan event range but not registered";
        LOG(DEBUG) << "(Shifted) Trigger ID " << triggerN << " within Corryvreckan event range but not registered";
    } else {
        // Redefine EUDAQ2 event begin and end to trigger timestamp (both were zero):
        evt->SetTimeBegin(static_cast<uint64_t>(clipboard->getEvent()->getTriggerTime(triggerN) * 1000));
        evt->SetTimeEnd(static_cast<uint64_t>(clipboard->getEvent()->getTriggerTime(triggerN) * 1000));
        LOG(DEBUG) << "Shifted Trigger ID " << triggerN << " found in Corryvreckan event";
    }
    return trigger_position;
}

int AnalysisItkStripEfficiency::get_next_tdc(const std::shared_ptr<Clipboard>& clipboard,
                                             const eudaq::FileReaderUP& filereader) {
    LOG(DEBUG) << "Get next event.";
    auto stdevt = eudaq::StandardEvent::MakeShared();
    while(1) {
        auto evt = filereader->GetNextEvent();
        //    LOG(INFO) << "Event Number being processed: " << evt->GetEventN();
        if(!evt) {
            LOG(DEBUG) << "Reached end-of-file.";
            throw EndOfFile();
        }
        LOG(DEBUG) << "Converting to StdEvent";
        if(!eudaq::StdEventConverter::Convert(evt, stdevt, eudaq_config_)) {
            return -1;
        }
        auto current_position = is_within_event(clipboard, stdevt);

        if(current_position == Event::Position::BEFORE) {
            LOG(DEBUG) << "Before current event, searching on";
            continue;
        }
        break;
    }
    LOG(DEBUG) << "Getting Tag " << ttc_tag_;
    if(!stdevt->HasTag(ttc_tag_)) {
        LOG(ERROR) << "Event tag (" << ttc_tag_ << ") is not available in the event. Check your data.";
        return -1;
    }

    int ttc_tag_content = stdevt->GetTag(ttc_tag_, -1);
    int delay = get_highest_bit_set(ttc_tag_content);
    LOG(DEBUG) << "Tag content:" << ttc_tag_content << "; Delay value: " << delay;
    return delay;
}

StatusCode AnalysisItkStripEfficiency::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Get the telescope tracks from the clipboard
    auto tracks = clipboard->getData<Track>();

    // auto pitch_x = m_detector->getPitch().X();
    // auto pitch_y = m_detector->getPitch().Y();
    pitch_x = m_detector->getPitch().X();
    pitch_y = m_detector->getPitch().Y();

    int delay = 0; // this is your new favourite value!

    try {
        delay = get_next_tdc(clipboard, readerTTC_);
    } catch(EndOfFile&) {
        return StatusCode::EndRun;
    }

    // Loop over all tracks
    for(auto& track : tracks) {
        n_track++;
        bool has_associated_cluster = false;
        bool is_within_roi = true;
        bool is_in_delay_window = true;
        LOG(DEBUG) << "Looking at next track";

        // Jonas making residual vs time plot before filtering events away
        // Get the event:
        // auto event = clipboard->getEvent();

        // Get the DUT clusters from the clipboard, that are assigned to the track
        // auto associated_clusters = track->getAssociatedClusters(m_detector->getName());

        // for(auto& associated_cluster : associated_clusters) {
        //     // Get the track intercept with the detector
        //     auto trackIntercept = m_detector->getIntercept(track.get());

        //     // Calculate the local residuals
        //     auto residual2D = m_detector->Residual(trackIntercept, associated_cluster);

        //      hresi_x_track_time->Fill(track->timestamp(),
        //         static_cast<double>(Units::convert(residual2D.X(), "mrad")));

        //      hresi_x_event_end->Fill(static_cast<double>(Units::convert(event->end(),"s")),
        //         static_cast<double>(Units::convert(residual2D.X(), "mrad")));

        //      hresi_x_event_Start->Fill(static_cast<double>(Units::convert(event->end(),"s")),
        //         static_cast<double>(Units::convert(residual2D.X(), "mrad")));

        //  }
        // // Cut on the chi2/ndof
        if(track->getChi2ndof() > m_chi2ndofCut) {
            LOG(DEBUG) << " - track discarded due to Chi2/ndof";
            n_chi2++;
            continue;
        }

        // Check if it intercepts the DUT
        auto globalIntercept = m_detector->getIntercept(track.get());
        auto localIntercept = m_detector->globalToLocal(globalIntercept);

        LOG(TRACE) << " Checking if track is outside DUT area";
        if(!m_detector->hasIntercept(track.get(), m_perimeter_exclude)) {
            LOG(DEBUG) << " - track outside DUT area: " << localIntercept;
            n_dut++;
            continue;
        }

        // Check that track is within region of interest using winding number algorithm
        LOG(TRACE) << " Checking if track is outside ROI";
        if(!m_detector->isWithinROI(track.get())) {
            LOG(DEBUG) << " - track outside ROI";
            n_roi++;
            is_within_roi = false;
            // here we don't continue because only some particular histograms shall be effected
        }

        // Check that it doesn't go through/near a masked pixel
        LOG(TRACE) << " Checking if track is close to masked pixel";
        if(m_detector->hitMasked(track.get(), m_maskedPixelDistanceCut)) {
            n_masked++;
            LOG(DEBUG) << " - track close to masked pixel";
            continue;
        }

        // Check that track is within delay window
        // delay < minimum  || delay >= maximum then False
        if((delay < m_delay_cuts[1]) || (delay >= m_delay_cuts[0])) {
            LOG(DEBUG) << "- track outside the delay window";
            LOG(DEBUG) << "- track delay: " << delay;
            LOG(DEBUG) << "- delay window from: " << m_delay_cuts[1] << " - " << m_delay_cuts[0];
            is_in_delay_window = false;
            n_timingWindow++;
        }

        // Get the event:
        auto event = clipboard->getEvent();

        // Discard tracks which are very close to the frame edges
        if(fabs(track->timestamp() - event->end()) < m_timeCutFrameEdge) {
            // Late edge - eventEnd points to the end of the frame`
            LOG(INFO) << " - track close to end of readout frame: "
                      << Units::display(fabs(track->timestamp() - event->end()), {"us", "ns"}) << " at "
                      << Units::display(track->timestamp(), {"us"});
            n_frameedge++;
            continue;
        } else if(fabs(track->timestamp() - event->start()) < m_timeCutFrameEdge) {
            // Early edge - eventStart points to the beginning of the frame
            LOG(INFO) << " - track close to start of readout frame: "
                      << Units::display(fabs(track->timestamp() - event->start()), {"us", "ns"}) << " at "
                      << Units::display(track->timestamp(), {"us"});
            n_frameedge++;
            continue;
        }

        // check if track has an associated cluster on required detector(s):
        auto foundRequiredAssocCluster = [this](Track* t) {
            for(auto& requireAssocCluster : require_associated_cluster_on_) {
                if(!requireAssocCluster.empty() && t->getAssociatedClusters(requireAssocCluster).size() == 0) {
                    LOG(DEBUG) << "No associated cluster from required detector " << requireAssocCluster << " on the track.";
                    return false;
                }
            }
            return true;
        };
        if(!foundRequiredAssocCluster(track.get())) {
            n_requirecluster++;
            continue;
        }

        // Count this as reference track:
        if(is_in_delay_window) {
            total_tracks++;
        }

        // Calculate in-pixel position of track in microns
        auto inpixel = m_detector->inPixel(localIntercept);
        auto xmod = inpixel.X();
        auto ymod = inpixel.Y();
        auto xmod_um = xmod * 1000.; // mm->um (for plotting)
        auto ymod_um = ymod * 1000.; // mm->um (for plotting)

        bool isWithinInPixelROI =
            (pitch_x - abs(xmod * 2) > m_inpixelEdgeCut.x()) && (pitch_y - abs(ymod * 2) > m_inpixelEdgeCut.y());

        // Get the DUT clusters from the clipboard, that are assigned to the track
        auto associated_clusters = track->getAssociatedClusters(m_detector->getName());
        if(associated_clusters.size() > 0) {
            hPos_TrackLocal_Ass->Fill(localIntercept.x(), localIntercept.y());

            auto cluster = track->getClosestCluster(m_detector->getName());
            has_associated_cluster = true;
            if(is_in_delay_window) {
                matched_tracks++;
            }
            auto pixels = cluster->pixels();
            for(auto& pixel : pixels) {
                if((pixel->column() == static_cast<int>(m_detector->getColumn(localIntercept)) &&
                    pixel->row() == static_cast<int>(m_detector->getRow(localIntercept))) &&
                   isWithinInPixelROI) {
                    hPixelEfficiencyMatrix_TProfile->Fill(pixel->column(), pixel->row(), 1);
                    break; // There cannot be a second pixel within the cluster through which the track goes.
                }
            }

            auto clusterLocal = m_detector->globalToLocal(cluster->global());

            auto distance =
                ROOT::Math::XYZVector(localIntercept.x() - clusterLocal.x(), localIntercept.y() - clusterLocal.y(), 0);
            hDistanceCluster_track->Fill(distance.X(), distance.Y());
            hDistanceCluster->Fill(std::sqrt(distance.Mag2()));

            hGlobalEfficiencyMap_clustPos_TProfile->Fill(
                cluster->global().x(), cluster->global().y(), has_associated_cluster);
            hGlobalEfficiencyMap_clustPos->Fill(has_associated_cluster, cluster->global().x(), cluster->global().y());

            hChipEfficiencyMap_clustPos_TProfile->Fill(
                m_detector->getColumn(clusterLocal), m_detector->getRow(clusterLocal), has_associated_cluster);
            hChipEfficiencyMap_clustPos->Fill(
                has_associated_cluster, m_detector->getColumn(clusterLocal), m_detector->getRow(clusterLocal));
        } else {
            hPos_TrackLocal_No_Ass->Fill(localIntercept.x(), localIntercept.y());
        }

        if(!has_associated_cluster && isWithinInPixelROI) {
            hPixelEfficiencyMatrix_TProfile->Fill(
                m_detector->getColumn(localIntercept), m_detector->getRow(localIntercept), 0);
        }

        hGlobalEfficiencyMap_trackPos_TProfile->Fill(globalIntercept.X(), globalIntercept.Y(), has_associated_cluster);
        hGlobalEfficiencyMap_trackPos->Fill(has_associated_cluster, globalIntercept.X(), globalIntercept.Y());

        hChipEfficiencyMap_trackPos_TProfile->Fill(
            m_detector->getColumn(localIntercept), m_detector->getRow(localIntercept), has_associated_cluster);
        hChipEfficiencyMap_trackPos->Fill(
            has_associated_cluster, m_detector->getColumn(localIntercept), m_detector->getRow(localIntercept));

        // For pixels, only look at the ROI:

        /* COMPLETELY WRONG FOR RADIAL
        // track fitting provides sub-strip positional resolution.
        localIntercept/pitch gets strip col id in floating form
        localIntercept/pitch - floor(localIntercept/pitch) gives position within strip
         multiply factor *N gives position with N strips
        */
        auto x_remainder = ((localIntercept.x() / (m_detector->getPitch().X() * 2.0)) -
                            (floor(localIntercept.x() / (m_detector->getPitch().X() * 2.0)))) *
                           2.0;
        auto y_remainder = ((localIntercept.y() / (m_detector->getPitch().Y() * 2.0)) -
                            (floor(localIntercept.y() / (m_detector->getPitch().Y() * 2.0)))) *
                           2.0;
        auto polar_det = std::dynamic_pointer_cast<PolarDetector>(m_detector);
        if(polar_det != nullptr) {
            auto localInterceptPolar = polar_det->getPolarPosition(localIntercept);
            x_remainder = ((localInterceptPolar.x() / (m_detector->getPitch().X() * 2.0)) -
                           (floor(localInterceptPolar.x() / (m_detector->getPitch().X() * 2.0)))) *
                          2.0;
            y_remainder = ((localInterceptPolar.y() / (m_detector->getPitch().Y() * 2.0)) -
                           (floor(localInterceptPolar.y() / (m_detector->getPitch().Y() * 2.0)))) *
                          2.0;
        }

        if(is_within_roi) {
            LOG(DEBUG) << "is_within_roi True, filling eTimingEfficiency";
            eTimingEfficiency->Fill(has_associated_cluster, delay);

            if(delay < 5 || delay > 28) {
                LOG(INFO) << "eTimingEfficiency was filled with " << has_associated_cluster << "  " << delay;
            }

            // LOG(DEBUG) << "eTimingEfficiency pointer print " << eTimingEfficiency;

            LOG(DEBUG) << "(has_associated_cluster, delay, x_remainder)" << has_associated_cluster << " , " << delay << ","
                       << x_remainder;

            hStripEfficiencyOddEven_TProfile->Fill(int(m_detector->getColumn(localIntercept)) % 2,
                                                   m_detector->getPitch().Y() * (y_remainder - floor(y_remainder)),
                                                   has_associated_cluster);
            // std::cout << "xr: " << x_remainder << " ylocal: " << localIntercept.y() << " row: " <<
            // m_detector->getRow(localIntercept) << std::endl;
            if(is_in_delay_window) {
                LOG(DEBUG) << " is_within_roi & is_in_delay_window, filling eTotalEfficiency , " << has_associated_cluster;
                eInPixelEfficiency->Fill(has_associated_cluster, x_remainder);
                eTotalEfficiency->Fill(has_associated_cluster, 0); // use 0th bin for total efficiency
                hPixelEfficiencyMap_trackPos_TProfile->Fill(xmod_um, ymod_um, has_associated_cluster);
                hPixelEfficiencyMap_trackPos->Fill(has_associated_cluster, xmod_um, ymod_um);
                efficiencyColumns->Fill(has_associated_cluster, m_detector->getColumn(localIntercept));
                efficiencyRows->Fill(has_associated_cluster, m_detector->getRow(localIntercept));
                efficiencyVsTime->Fill(has_associated_cluster, track->timestamp() / 1e9); // convert nanoseconds to seconds
                LOG(DEBUG) << "efficiencyVsTime filled with: " << has_associated_cluster << " " << track->timestamp() / 1e9;

                if(isWithinInPixelROI) {
                    LOG(DEBUG) << "isWithinInPixelROI true, filling eTotalEfficiency_inPixelROI" << xmod_um << " " << ymod_um
                               << " " << has_associated_cluster;
                    hPixelEfficiencyMap_inPixelROI_trackPos_TProfile->Fill(xmod_um, ymod_um, has_associated_cluster);
                    eTotalEfficiency_inPixelROI->Fill(has_associated_cluster, 0); // use 0th bin for total efficiency
                }
            }
        }

        auto intercept_col = static_cast<size_t>(m_detector->getColumn(localIntercept));
        auto intercept_row = static_cast<size_t>(m_detector->getRow(localIntercept));

        if(has_associated_cluster) {
            hTimeDiffPrevTrack_assocCluster->Fill(
                static_cast<double>(Units::convert(track->timestamp() - last_track_timestamp, "us")));
            hRowDiffPrevTrack_assocCluster->Fill(m_detector->getRow(localIntercept) - last_track_row);
            hColDiffPrevTrack_assocCluster->Fill(m_detector->getColumn(localIntercept) - last_track_col);
            hPosDiffPrevTrack_assocCluster->Fill(m_detector->getColumn(localIntercept) - last_track_col,
                                                 m_detector->getRow(localIntercept) - last_track_row);
            if((prev_hit_ts.at(intercept_col)).at(intercept_row) != 0) {
                hTrackTimeToPrevHit_matched->Fill(static_cast<double>(
                    Units::convert(track->timestamp() - prev_hit_ts.at(intercept_col).at(intercept_row), "us")));
            }
        } else {
            hGlobalEfficiencyMap_clustPos_TProfile->Fill(globalIntercept.X(), globalIntercept.Y(), has_associated_cluster);
            hGlobalEfficiencyMap_clustPos->Fill(has_associated_cluster, globalIntercept.X(), globalIntercept.Y());

            hChipEfficiencyMap_clustPos_TProfile->Fill(
                m_detector->getColumn(localIntercept), m_detector->getRow(localIntercept), has_associated_cluster);
            hChipEfficiencyMap_clustPos->Fill(
                has_associated_cluster, m_detector->getColumn(localIntercept), m_detector->getRow(localIntercept));

            hTimeDiffPrevTrack_noAssocCluster->Fill(
                static_cast<double>(Units::convert(track->timestamp() - last_track_timestamp, "us")));
            hRowDiffPrevTrack_noAssocCluster->Fill(m_detector->getRow(localIntercept) - last_track_row);
            hColDiffPrevTrack_noAssocCluster->Fill(m_detector->getColumn(localIntercept) - last_track_col);
            hPosDiffPrevTrack_noAssocCluster->Fill(m_detector->getColumn(localIntercept) - last_track_col,
                                                   m_detector->getRow(localIntercept) - last_track_row);
            if((prev_hit_ts.at(intercept_col)).at(intercept_row) != 0) {
                LOG(DEBUG) << "Found a time difference of "
                           << Units::display(track->timestamp() - prev_hit_ts.at(intercept_col).at(intercept_row), "us");
                hTrackTimeToPrevHit_notmatched->Fill(static_cast<double>(
                    Units::convert(track->timestamp() - prev_hit_ts.at(intercept_col).at(intercept_row), "us")));
            }
        }
        last_track_timestamp = track->timestamp();
        last_track_col = m_detector->getColumn(localIntercept);
        last_track_row = m_detector->getRow(localIntercept);
    } // end loop over tracks

    // Before going to the next event, loop over all pixels (all hits incl. noise)
    // and fill matrix with timestamps of previous pixels.
    auto pixels = clipboard->getData<Pixel>(m_detector->getName());
    if(pixels.empty()) {
        LOG(DEBUG) << "Detector " << m_detector->getName() << " does not have any pixels on the clipboard";
    }

    for(auto& pixel : pixels) {
        if(pixel->column() < 0 || pixel->column() > m_detector->nPixels().X() || pixel->row() < 0 ||
           pixel->row() > m_detector->nPixels().Y()) {
            continue;
        }
        prev_hit_ts.at(static_cast<size_t>(pixel->column())).at(static_cast<size_t>(pixel->row())) = pixel->timestamp();
    }

    return StatusCode::Success;
}

void AnalysisItkStripEfficiency::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
    // Track selection flow:
    LOG(STATUS) << "Track selection flow:       " << n_track << std::endl
                << "* rejected by chi2          -" << n_chi2 << std::endl
                << "* track outside ROI         -" << n_roi << std::endl
                << "* track outside DUT         -" << n_dut << std::endl
                << "* track outside TOA Winwow  -" << n_timingWindow << std::endl
                << "* track close to masked px  -" << n_masked << std::endl
                << "* track close to frame edge -" << n_frameedge << std::endl
                << "* track without an associated cluster on required detector - " << n_requirecluster << std::endl
                << "Accepted tracks:            " << total_tracks;

    double totalEff = 100 * static_cast<double>(matched_tracks) / (total_tracks > 0 ? total_tracks : 1);
    double lowerEffError = totalEff - 100 * (TEfficiency::ClopperPearson(total_tracks, matched_tracks, 0.683, false));
    double upperEffError = 100 * (TEfficiency::ClopperPearson(total_tracks, matched_tracks, 0.683, true)) - totalEff;
    LOG(STATUS) << "Total efficiency of detector " << m_detector->getName() << ": " << totalEff << "(+" << upperEffError
                << " -" << lowerEffError << ")%, measured with " << matched_tracks << "/" << total_tracks
                << " matched/total tracks";

    for(int icol = 1; icol < m_detector->nPixels().X() + 1; icol++) {
        for(int irow = 1; irow < m_detector->nPixels().Y() + 1; irow++) {
            // calculate total efficiency: (just to double check the other calculation)
            const int bin = hChipEfficiencyMap_trackPos->GetGlobalBin(icol, irow);
            double eff = hChipEfficiencyMap_trackPos->GetEfficiency(bin);
            if(eff > 0) {
                LOG(TRACE) << "col/row = " << icol << "/" << irow << ", binContent = " << eff;
                hPixelEfficiency->Fill(eff);
            }
            eff = hPixelEfficiencyMatrix_TProfile->GetBinContent(bin);
            if(eff > 0) {
                LOG(TRACE) << "col/row = " << icol << "/" << irow << ", binContent = " << eff;
                hPixelEfficiencyMatrix->Fill(eff);
            }
        }
    }
}
