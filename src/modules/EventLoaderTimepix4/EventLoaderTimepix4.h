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
        enum uint8_t {
            pixel_data = 0x00,

            ctrl_heartbeat = 0xE0,
            shutter_rise = 0xE1,
            shutter_fall = 0xE2,
            t0_sync = 0xE3,
            signal_rise = 0xE4,
            signal_fall = 0xE5,
            ctrl_data_test = 0xEA,

            frame_start = 0xF0,
            frame_end = 0xF1,
            segment_start = 0xF2,
            segment_end = 0xF3,

            header_invalid = 0xFF
        };

        struct pixelData {
            uint64_t fullTot;
            uint64_t fullToa;
            uint64_t hb;
            bool t0;
            bool isDigital;
        };

        struct heartbeatData {
            uint64_t time;
            uint64_t bufferID;
        };

        std::shared_ptr<Detector> m_detector;

        // ROOT graphs
        TH2F* hHitMap;
        TH1F* hRawToT;
        TH1F* hRawFullToT;
        TH1F* hToT;
        TH1F* hRawToA;
        TH1F* hRawExtendedToA;
        TH1F* hRawFullToA;
        TH1F* hHitTime;

        bool decodeNextWord();
        bool decodePacket(uint64_t packet);
        void fillBuffer();
        bool loadData(const std::shared_ptr<Clipboard>& clipboard, PixelVector&);

        // configuration parameters:
        std::string m_inputDirectory;

        std::string calibrationPath;
        std::string threshold;

        std::vector<std::vector<float>> vtot;
        std::vector<std::vector<float>> vtoa;

        heartbeatData m_hbData;
        uint16_t m_hbIndex = 0;
        std::vector<heartbeatData> m_hbDataBuffer;

        std::vector<pixelData> m_pDataBuffer;




        // pixel packet variables
        uint64_t m_addr;
        uint64_t m_pileup;
        uint64_t m_tot;
        uint64_t m_ftoa_fall;
        uint64_t m_ftoa_rise;
        uint64_t m_uftoa_start;
        uint64_t m_uftoa_stop;
        uint64_t m_ext_toa;
        uint16_t m_toa;
        uint64_t m_pixel;
        uint64_t m_sPixel;
        uint64_t m_sPGroup;
        uint64_t m_fullTot;
        uint64_t m_fullToa;
        uint64_t m_oldbeat;


        uint64_t m_packetTime[2] = {0};

        std::tuple<uint32_t, uint32_t> m_colrow;

        // conversion factor from:
        // heartbeat tdc (25ns)
        // fine tdc (~1.56ns | 1/640 MHz^-1)
        // ultrafine tdc (~195 ps | 1/(8*640) MHz^-1)

        // location of the digital pixels
        std::tuple<uint32_t, uint32_t> m_digColRow[8] = {
            {0, 0}, {4, 1}, {441, 2}, {445, 3}, {2, 508}, {6, 509}, {443, 510}, {447, 511}};

        std::streampos m_stream_pos[2] = {0, 0};
        uint m_fIndex = 0;

        // Member variables
        std::vector<std::unique_ptr<std::ifstream>> m_files;
        std::vector<std::unique_ptr<std::ifstream>>::iterator m_file_iterator;


        unsigned long long int m_syncTime {0};
        uint64_t m_unsynced[2] = {1, 1};
        bool m_clearedHeader {false};
        long long int m_syncTimeTDC {0};
        int m_TDCoverflowCounter {0};
        int m_prevTriggerNumber {0};
        int m_triggerOverflowCounter {0};
        bool eof_reached {false};

        size_t m_buffer_depth;
        std::vector<uint64_t> m_dataBuffer;
        uint64_t m_dataPacket;
        long long int m_currentEvent;


        template <typename T> struct CompareTimeGreater {
            bool operator()(const std::shared_ptr<T> a, const std::shared_ptr<T> b) {
                return a->timestamp() > b->timestamp();
            }
        };

        std::priority_queue<std::shared_ptr<Pixel>, PixelVector, CompareTimeGreater<Pixel>> sorted_pixels_;

        //======================================================================================================================
        // Unpack the header of the data packet (original, kepler)
        //======================================================================================================================
        std::array<unsigned, 5> decode_header(uint64_t packet) const {
            return {unsigned(0xF & (packet >> 60)),
                    unsigned(0x3 & (packet >> 58)),
                    unsigned(0x3FF & (packet >> 48)),
                    unsigned(0x1FFF & (packet >> 32)),
                    unsigned(0xFFFFFFFF & (packet >> 0))};
        }

        uint64_t fullTot(uint64_t ftoa_rise, uint64_t ftoa_fall, uint64_t uftoa_start, uint64_t uftoa_stop, uint64_t tot) {
            // Decode TOT. Units are period of 8*640MHz (195 ps)
            return ((tot << 7) + ((ftoa_rise - ftoa_fall) << 3) - (uftoa_start - uftoa_stop));
        }

        uint64_t fullToa(uint64_t toa, uint64_t uftoa_start, uint64_t uftoa_stop, uint64_t ftoa_rise) {
            // Decode TOA. Units are period of 8*640MHz (195 ps)
            return ((toa << 7) - (ftoa_rise << 3) + (uftoa_start - uftoa_stop));
        }

        uint64_t toa_clkdll_correction(uint64_t spgroup_addr = 0) {
            // Corrects latency delay due to DDLL clock distribution. Units are period of 40MHz (25ns)
            uint64_t clk_dll_step = 1 >> 5;
            return (15 - spgroup_addr) * clk_dll_step;
        }

        // decodes the row and column position from the address dat etc.
        std::tuple<uint32_t, uint32_t>
        decodeColRow(uint64_t pix, uint64_t sPix, uint64_t spixgrp, uint64_t header, bool top) { // taken from spidr4tools
            uint32_t col;
            uint32_t row;
            col = static_cast<uint32_t>(header << 1 | pix >> 2);
            row = static_cast<uint32_t>(spixgrp << 4 | sPix << 2 | (pix & 0x3));
            if(top) // top half counting is inverted
            {
                col = 448 - 1 - col;
                row = 512 - 1 - row;
            }

            return {col, row};
        }

        uint64_t getAddr(uint64_t packet) {
            return (packet >> 46) & 0x3ffff;
        } // address including pixel, super pixel and super pixel group values
        uint64_t getSuperPixelGroup(uint64_t packet) { return (packet >> 51) & 0xf; } // super pixel group address
        uint64_t getSuperPixel(uint64_t packet) { return (packet >> 49) & 0x3; }      // super pixel address
        uint64_t getPixel(uint64_t packet) { return (packet >> 46) & 0x7; }           // pixel address

        uint16_t getToA(uint64_t packet) {
            return (packet >> 30) & 0xffff;
        } // Time of Arrival (ToA) | units of 25 ns (1/(40 MHz))
        uint64_t getFToARise(uint64_t packet) {
            return (packet >> 17) & 0x1f;
        } // fine ToA rising edge | units of ~1.56 ns (1/(640 MHz)
        uint64_t getFToAFall(uint64_t packet) {
            return (packet >> 12) & 0x1f;
        } // fine ToA falling edge | units of ~1.56 ns (1/(640 MHz)
        uint64_t getToT(uint64_t packet) {
            return (packet >> 1) & 0x7ff;
        } // Time over Threshold | units of 25 ns  (1/(40 MHz))
        uint64_t getPileUp(uint64_t packet) {
            return packet & 0x1;
        } // bit to register whether another hit was coming in while pixel was busy

        bool compareTupleEq(std::tuple<uint32_t, uint32_t> tuple1, std::tuple<uint32_t, uint32_t> tuple2) {
            if(std::get<0>(tuple1) == std::get<0>(tuple2) && std::get<1>(tuple1) == std::get<1>(tuple2)) {
                return true;
            } else {
                return false;
            }
        }

        std::tuple<uint, std::vector<std::unique_ptr<std::ifstream>>::iterator>
        switchHalf(uint fIndex, std::vector<std::unique_ptr<std::ifstream>>::iterator fIterator) {
            if(fIndex) {
                fIterator--;
                fIndex = 0;
                return {fIndex, fIterator};
            } else {
                fIterator++;
                fIndex = 1;
                return {fIndex, fIterator};
            }
        }

        uint64_t extendToa(uint64_t toa, uint64_t heartbeat, uint64_t tot) {
            // extending toa by heartbeat counter
            uint64_t extToa = toa | (heartbeat & 0xFFFFFFFFFFFF0000);

            // toa vs heartbeat for latency correction
            if(extToa + 0x8000 < heartbeat)
                toa += 0x10000;
            else if(extToa > heartbeat + 0x8000 && toa >= 0x10000)
                toa -= 0x10000;

            if(!tot)
                extToa++;
            return extToa;
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
        // uftoa encoding change to actual values (original, TPX4TOOLS by kevin)
        // =======================================
        uint64_t uftoaBin[16] = {4, 5, 8, 6, 8, 8, 8, 7, 3, 8, 8, 8, 2, 8, 1, 0};
        uint64_t UftoaStart(uint64_t value) { return uftoaBin[value >> 26 & 0x000F]; }
        uint64_t UftoaStop(uint64_t value) { return uftoaBin[value >> 22 & 0x000F]; }
    };
} // namespace corryvreckan
#endif // Timepix4EVENTLOADER_H
