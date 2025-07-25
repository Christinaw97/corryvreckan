/**
 * @file
 * @brief Definition of module ClusteringSpatial
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef ClusteringSpatial_H
#define ClusteringSpatial_H 1

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class ClusteringSpatial : public Module {

    public:
        // Constructors and destructors
        ClusteringSpatial(Configuration& config, std::shared_ptr<Detector> detector);
        ~ClusteringSpatial() {}

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        std::shared_ptr<Detector> m_detector;

        void calculateClusterCentre(Cluster*);

        // Cluster histograms
        TH1F* clusterSize;
        TH1F* clusterSeedCharge;
        TH1F* clusterWidthRow;
        TH1F* clusterWidthColumn;
        TH1F* clusterCharge;
        TH1F* clusterMultiplicity;
        TH2F* clusterPositionGlobal;
        TH2F* clusterPositionLocal;
        TH1F* clusterTimes;
        TH1F* clusterUncertaintyX;
        TH1F* clusterUncertaintyY;

        bool useTriggerTimestamp;
        bool chargeWeighting;
        bool rejectByROI;
        int neighbor_radius_row_;
        int neighbor_radius_col_;
    };
} // namespace corryvreckan
#endif // ClusteringSpatial_H
