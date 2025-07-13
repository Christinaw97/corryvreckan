/**
 * @file
 * @brief Definition of module TreeWriter
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
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

typedef struct 
{
  double    xSlope;
  double    ySlope;
  double    xIntercept;
  double    yIntercept;
  double    chi2;
  double    xResidBack;
  double    yResidBack;
  double    xErrDUT;
  double    yErrDUT;
  double    xErr04;
  double    yErr04;
  double    xErr05;
  double    yErr05;
  double    xErrPix0;
  double    yErrPix0;
  double    xResid04;
  double    yResid04;
  double    xResid05;
  double    yResid05;
  double    xResidPix0;
  double    yResidPix0;
  int       trigger;
  int       runNumber;
  int       nPlanes;
  int       numPixels;
  int       numBackPlanes;
  int       numTracks;
  int       numClustersPix;
  int       numClustersStripsOdd;
  int       numClustersStripsEven;
  int       numStripsWith2Clusters;
  long long timestamp;
  long long bco;
} ConvertedEvent;

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class TreeWriter : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detectors Vector of pointers to the detectors
         */
        TreeWriter(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors);

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

        // // Output data file to write
        // std::unique_ptr<TFile> output_file_;
        // std::string output_file_name_{};

        // // List of trees that are stored in data file
        // std::map<std::string, std::unique_ptr<TTree>> trees_;
        // std::unique_ptr<TTree> event_tree_;
        // // Event* event_{};

         std::shared_ptr<Detector> m_detector;




        // double trackChi2;
        // size_t trackNClusters;
        // size_t ntracks;
        // size_t trackNdof;
        // double yIntercept;
        // double xIntercept;
        // double xSlope;
        // double ySlope;

        ConvertedEvent theConvertedEvent_;

       
        std::map<std::string, Object*> m_objects;

        TFile* m_outputFile;
        TTree* m_outputTree{};

        // Config parameters
        std::string m_fileName;
        std::string m_treeName;
    };

} // namespace corryvreckan
