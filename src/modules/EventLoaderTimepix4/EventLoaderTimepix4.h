/**
 * @file
 * @brief Definition of module EventLoaderTimepix4
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef TIMEPIX4EVENTLOADER_H
#define TIMEPIX4EVENTLOADER_H 1

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <queue>
#include <stdio.h>
#include "core/module/Module.hpp"
#include "objects/Pixel.hpp"
#include "objects/SpidrSignal.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class EventLoaderTimepix4 : public Module {

    public:
        // Constructors and destructors
        EventLoaderTimepix4(Configuration& config, std::shared_ptr<Detector> detector);
        ~EventLoaderTimepix4() {}

        // Standard algorithm functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        std::shared_ptr<Detector> m_detector;
        long m_numOfEvents;

        // ROOT graphs
        TH2F* hHitMap;
        TH1F* pixelToT_beforecalibration;
        TH1F* pixelToT_aftercalibration;
        TH2F* pixelTOTParameterA;
        TH2F* pixelTOTParameterB;
        TH2F* pixelTOTParameterC;
        TH2F* pixelTOTParameterT;
        TH2F* pixelTOAParameterC;
        TH2F* pixelTOAParameterD;
        TH2F* pixelTOAParameterT;
        TH1F* timeshiftPlot;
        TH1F* hTriggerTime;


        bool decodeNextWord();
        bool decodePacket(uint64_t packet, uint64_t ratio_VCO_CKDLL, bool gray);
        void fillBuffer();
        bool loadData(const std::shared_ptr<Clipboard>& clipboard, PixelVector&, SpidrSignalVector&);

        // configuration parameters:
        std::string m_inputDirectory;

        bool applyCalibration;
        std::string calibrationPath;
        std::string threshold;

        std::vector<std::vector<float>> vtot;
        std::vector<std::vector<float>> vtoa;


        // pixel packet variables
        uint64_t m_addr;
        uint64_t m_pileup;
        uint64_t m_tot;
        uint64_t m_ftoa_fall;
        uint64_t m_ftoa_rise;
        uint64_t m_uftoa_start;
        uint64_t m_uftoa_stop;
        uint64_t m_toa;
        uint64_t m_pixel;
        uint64_t m_sPixel;
        uint64_t m_sPGroup;
        uint64_t m_fullTot;
        uint64_t m_fullToa;
        uint64_t m_heartbeat;
        uint64_t m_oldside;
        std::tuple<uint32_t, uint32_t> m_colrow;

        std::streampos m_stream_pos[2] = {0};
        uint m_file_index = 0;


        // Member variables
        std::vector<std::unique_ptr<std::ifstream>> m_files;
        std::vector<std::unique_ptr<std::ifstream>>::iterator m_file_iterator;

        bool eof_reached;
        size_t m_buffer_depth;
        std::vector<uint64_t> m_dataBuffer;
        uint64_t m_dataPacket;
        unsigned long long int m_syncTime;
        bool m_clearedHeader;
        long long int m_syncTimeTDC;
        int m_TDCoverflowCounter;



        long long int m_currentEvent;

        unsigned long long int m_prevTime;
        bool m_shutterOpen;
        int m_prevTriggerNumber;
        int m_triggerOverflowCounter;

        template <typename T> struct CompareTimeGreater {
            bool operator()(const std::shared_ptr<T> a, const std::shared_ptr<T> b) {
                return a->timestamp() > b->timestamp();
            }
        };

        std::priority_queue<std::shared_ptr<Pixel>, PixelVector, CompareTimeGreater<Pixel>> sorted_pixels_;
        std::priority_queue<std::shared_ptr<SpidrSignal>, SpidrSignalVector, CompareTimeGreater<SpidrSignal>>
            sorted_signals_;

        //======================================================================================================================
        // Unpack the header of the data packet (original, kepler)
        //======================================================================================================================
        std::array<unsigned, 5> decode_header(uint64_t packet) const {
          return {unsigned(0xF & (packet >> 60)), unsigned(0x3 & (packet >> 58)), unsigned(0x3FF & (packet >> 48)),
                  unsigned(0x1FFF & (packet >> 32)), unsigned(0xFFFFFFFF & (packet >> 0))};
        }

        uint64_t decode_tot(uint64_t ftoa_rise, uint64_t ftoa_fall, uint64_t uftoa_start, uint64_t uftoa_stop, uint64_t tot, uint64_t ratio_VCO_CKDLL=16) {
            //Decode TOT. Units are period of 40MHz (25ns)
            return (tot + (ftoa_rise - ftoa_fall)/ratio_VCO_CKDLL - (uftoa_start-uftoa_stop)/(8 * ratio_VCO_CKDLL));
        }

        uint64_t decode_toa(uint64_t toa, uint64_t uftoa_start, uint64_t uftoa_stop, uint64_t ftoa_rise, uint64_t ratio_VCO_CKDLL=16){
          //Decode TOA. Units are period of 40MHz (25ns)
          return (toa - ftoa_rise/ratio_VCO_CKDLL + (uftoa_start-uftoa_stop)/(8 * ratio_VCO_CKDLL));
        }

        uint64_t toa_clkdll_correction(uint64_t spgroup_addr=0){
          //Corrects latency delay due to DDLL clock distribution. Units are period of 40MHz (25ns)
          uint64_t clk_dll_step=1/32;
          return (15-spgroup_addr)*clk_dll_step;
        }

        std::tuple<uint32_t, uint32_t> decodeColRow(uint64_t pix, uint64_t sPix, uint64_t spixgrp, uint64_t header, bool top){ // taken from spidr4tools
            uint32_t col;
            uint32_t row;
            col = header << 1 | pix >> 2;
            row = spixgrp << 4 | sPix << 2 | (pix & 0x3);
            if(top) // top half counting is inverted
            {
                col  = 448 - 1 - col;
                row  = 512 - 1 - row;
            }

            return {col,row};
        }

        inline uint16_t GrayToBin(uint16_t val) // taken from spidr4tools
        {
            val ^= val >> 8;
            val ^= val >> 4;
            val ^= val >> 2;
            val ^= val >> 1;

            return val;
        }

        //========================================
        // uftoa encoding change to acutal values (original, TPX4TOOLS by kevin)
        // =======================================
        uint64_t uftoaBin[16] = {4, 5, 8, 6, 8, 8, 8, 7, 3, 8, 8, 8, 2, 8, 1, 0};
        uint64_t UftoaStart(uint64_t value){ return uftoaBin[value >> 26 & 0x000F]; }
        uint64_t UftoaStop(uint64_t value){ return uftoaBin[value >> 22 & 0x000F]; }



    };
} // namespace corryvreckan
#endif // Timepix4EVENTLOADER_H
