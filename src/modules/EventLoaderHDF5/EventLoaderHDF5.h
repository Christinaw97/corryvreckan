/**
 * @file
 * @brief Definition of module EventLoaderHDF5
 *
 * @copyright Copyright (c) 2023-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef EventLoaderHDF5_H
#define EventLoaderHDF5_H 1

#include <H5Cpp.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>

#include <queue>
#include "core/module/Module.hpp"
#include "objects/Pixel.hpp"

namespace corryvreckan {

    class EventLoaderHDF5 : public Module {

    public:
        // Constructors and destructors
        EventLoaderHDF5(Configuration& config, std::shared_ptr<Detector> detector);
        ~EventLoaderHDF5() {}

        // Standard algorithm functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        struct Hit {
            int column;
            int row;
            int raw;
            double charge;
            double timestamp;
            uint32_t trigger_number;
        };

        using HitVector = std::vector<std::shared_ptr<Hit>>;

        H5::CompType h5_datatype{sizeof(Hit)};

        std::string m_fileName;
        std::string m_datasetName;
        double m_eventLength;
        std::shared_ptr<Detector> m_detector;
        hsize_t m_bufferDepth;
        bool m_sync_by_trigger;
        double m_timestampShift;
        uint32_t m_triggerShift;

        H5::DataSet m_dataset;
        H5::H5File m_file;
        H5::DataSpace f_dataspace;
        hsize_t f_total_records;
        hsize_t m_start_record;

        // Plots
        TH2F* hHitMap;
        TProfile2D* hTotMap;
        TH1F* hPixelToT;
        TH1F* hPixelCharge;
        TH1F* hPixelMultiplicity;
        TH1D* hClipboardEventStart;
        TH1D* hClipboardEventStart_long;
        TH1D* hClipboardEventEnd;
        TH1D* hClipboardEventDuration;

        // Additional helper function
        std::vector<Hit> readChunk();
        bool loadData(const std::shared_ptr<Clipboard>& clipboard, PixelVector&);
        void fillBuffer();
        Event::Position getPosition(const std::shared_ptr<Event>& event, const std::shared_ptr<Hit>& hit) const;

        std::priority_queue<std::shared_ptr<Hit>, HitVector, std::function<bool(std::shared_ptr<Hit>, std::shared_ptr<Hit>)>>
            m_buffer;
    };

} // namespace corryvreckan

#endif // EventLoaderHDF5_H
