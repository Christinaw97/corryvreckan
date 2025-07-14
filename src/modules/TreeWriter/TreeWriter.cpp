/**
 * @file
 * @brief Implementation of module TreeWriter
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "TreeWriter.h"

using namespace corryvreckan;

TreeWriter::TreeWriter(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)) {}

void TreeWriter::initialize() {    
    // // Initialise member variables
    LOG(DEBUG) << "Initialised TreeWriterDUT";

    config_.setDefault<std::string>("file_name", "outputTuples.root");

    m_fileName = config_.get<std::string>("file_name");

    // Create output file and directories
    auto path = createOutputFile(m_fileName, "root");
    m_outputFile = new TFile(path.c_str(), "RECREATE");
    LOG(DEBUG) << "Made and moved to output file: " << path;
    gDirectory->Delete("tree;*");

    m_outputTree = new TTree("CMSTiming", "The reconstructed telescope tracks");
    m_outputTree->Branch("event", &theConvertedEvent_, "xSlope/D:ySlope/D:xIntercept/D:yIntercept/D:chi2/D:xResidBack/D:yResidBack/D:xErrDUT/D:yErrDUT/D:xErr04/D:yErr04/D:xErr05/D:yErr05/D:xErrPix0/D:yErrPix0/D:xResid04/D:yResid04/D:xResid05/D:yResid05/D:xResidPix0/D:yResidPix0/D:trigger/I:runNumber/I:nPlanes/I:numPixels/I:numBackPlanes/I:numTracks/I:numClustersPix/I:numClustersStripsOdd/I:numClustersStripsEven/I:numStripsWith2Clusters:timestamp/L:bco/L");       



    LOG(DEBUG) << "Created tree: " << m_treeName;

    m_eventNumber = 0;

    // // Create the output branches
    // m_outputTree->Branch("i_evt", &m_eventNumber);
    // m_outputTree->Branch("ndof", &trackNdof);
    // m_outputTree->Branch("chi2", &trackChi2); //reduced chi2
    // m_outputTree->Branch("nClusters", &trackNClusters);
    // m_outputTree->Branch("ntracks", &ntracks);

    // m_outputTree->Branch("xIntercept", &xIntercept);
    // m_outputTree->Branch("yIntercept", &yIntercept);
    // m_outputTree->Branch("xSlope", &xSlope);
    // m_outputTree->Branch("ySlope", &ySlope);

  
}

StatusCode TreeWriter::run(const std::shared_ptr<Clipboard>& clipboard) {

    // //initial values
    theConvertedEvent_.chi2 = 999.;
    theConvertedEvent_.numClustersPix = 0.;
    theConvertedEvent_.xIntercept = -999.;
    theConvertedEvent_.yIntercept = -999.;
    theConvertedEvent_.xSlope = -999.;
    theConvertedEvent_.ySlope = -999.;
    theConvertedEvent_.trigger = m_eventNumber;

    // // Increment event counter
    m_eventNumber++;

    // // Getting tracks from the clipboard
    auto tracks = clipboard->getData<Track>();
    if (tracks.size() == 0) return StatusCode::Success;

    // // // Iterate through tracks found
    double chi2min = 999.;
    double chi2_temp = 999.;
    std::shared_ptr<corryvreckan::Track> minChi2Track = nullptr; // Pointer to the track with max chi2
    theConvertedEvent_.numTracks = tracks.size();
    for(auto& track : tracks) {
        chi2_temp = track->getChi2ndof();
        if (chi2_temp < chi2min) {
            chi2min = chi2_temp;
            minChi2Track = track;
        }   
    }
    if (minChi2Track) {
        theConvertedEvent_.chi2 = minChi2Track->getChi2ndof();
        theConvertedEvent_.numClustersPix = minChi2Track->getNClusters();
        ROOT::Math::XYZPoint intercept = minChi2Track->getIntercept(0.0);
        theConvertedEvent_.xIntercept = intercept.x()*1e3; //convert from mm to um
        theConvertedEvent_.yIntercept = intercept.y()*1e3; //convert from mm to um
        ROOT::Math::XYZVector slope = minChi2Track->getDirection(0.0);
        theConvertedEvent_.xSlope = slope.x();
        theConvertedEvent_.ySlope = slope.y();
    }
    

    
    m_outputTree->Fill();
    // Return value telling analysis to keep running
    return StatusCode::Success;
}

void TreeWriter::finalize(const std::shared_ptr<ReadonlyClipboard>&) { 
    LOG(DEBUG) << "Analysed " << m_eventNumber << " events"; 

    LOG(STATUS) << m_eventNumber << " events written to file " << m_fileName;

    // Writing out outputfile
    m_outputFile->Write();
    delete(m_outputFile);


}
