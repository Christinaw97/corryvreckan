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
#include <filesystem>
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
        enum headerIdentifier : uint8_t {
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
        TH1F* hFToARise;
        TH1F* hFToAFall;
        TH1F* hUfToAStop;
        TH1F* hUfToAStart;
        TH1F* hRawFullToA;
        TH1F* hHitTime;
        TH1F* hPileUp;

        // configuration parameters:
        std::string m_inputDirectory;
        std::filesystem::path m_inputPath;

        // data information
        heartbeatData m_hbData;
        uint16_t m_hbIndex = 0;
        std::vector<heartbeatData> m_hbDataBuffer;
        std::vector<pixelData> m_pDataBuffer;
        size_t m_buffer_depth;
        std::vector<uint64_t> m_dataBuffer;
        uint64_t m_dataPacket;
        long long int m_currentEvent;

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
        uint64_t m_packetTime[2] = {0, 0};

        // decoded column and row value
        std::tuple<uint32_t, uint32_t> m_colrow;

        // conversion factor from:
        // heartbeat tdc (25ns)
        // fine tdc (~1.56ns | 1/640 MHz^-1)
        // ultrafine tdc (~195 ps | 1/(8*640) MHz^-1)

        // location of the digital pixels in the matrix for filtering
        std::tuple<uint32_t, uint32_t> m_digColRow[8] = {
            {0, 0}, {4, 1}, {441, 2}, {445, 3}, {2, 508}, {6, 509}, {443, 510}, {447, 511}};

        uint m_fIndex = 0;

        // File Member variables
        std::vector<std::unique_ptr<std::ifstream>> m_files;
        std::vector<std::unique_ptr<std::ifstream>>::iterator m_file_iterator;

        // initialized variables for synchronization, header clearing etc.
        uint64_t m_unsynced[2] = {1, 1};
        bool eof_reached{false};

        //===============================================================
        // Begin of functions that are used within the Module
        //===============================================================
        template <typename T> struct CompareTimeGreater {
            bool operator()(const std::shared_ptr<T> a, const std::shared_ptr<T> b) {
                return a->timestamp() > b->timestamp();
            }
        };

        std::priority_queue<std::shared_ptr<Pixel>, PixelVector, CompareTimeGreater<Pixel>> sorted_pixels_;

        // decodes the next word, consisting of a header which tells me how many following datapackets are part of
        // (pixel/heartbeat data)
        bool decodeNextWord();

        // decodes the 64 bit dat apacket
        bool decodePacket(uint64_t packet);

        void fillBuffer();

        bool loadData(const std::shared_ptr<Clipboard>& clipboard, PixelVector&);

        // Unpack the header of the data packet. Taken from kepler
        std::array<unsigned, 5> decode_header(uint64_t packet);

        // decodes the row and column position from the address dat etc. Taken from spidr4tools
        std::tuple<uint32_t, uint32_t>
        decodeColRow(uint64_t pix, uint64_t sPix, uint64_t spixgrp, uint64_t header, bool top);

        // compares whether the tuple is part of the digital pixel array
        bool compareTupleEq(std::tuple<uint32_t, uint32_t> tuple1, std::tuple<uint32_t, uint32_t> tuple2);

        // switches the file iterator from one to the next
        std::tuple<uint, std::vector<std::unique_ptr<std::ifstream>>::iterator>
        switchHalf(uint fIndex, std::vector<std::unique_ptr<std::ifstream>>::iterator fIterator);

        // extension of the 16 bit ToA using the 64 bit heartbeat counter
        uint64_t extendToa(uint64_t toa, uint64_t heartbeat, uint64_t tot);

        // converts gray encoded bits to binary
        uint16_t GrayToBin(uint16_t val);

        // Decode TOT. Units are period of 8*640MHz (195 ps)
        uint64_t fullTot(uint64_t ftoa_rise, uint64_t ftoa_fall, uint64_t uftoa_start, uint64_t uftoa_stop, uint64_t tot) {
            return ((tot << 7) + ((ftoa_rise - ftoa_fall) << 3) - (uftoa_start - uftoa_stop));
        }

        // Decode TOA. Units are period of 8*640MHz (195 ps)
        uint64_t fullToa(uint64_t toa, uint64_t uftoa_start, uint64_t uftoa_stop, uint64_t ftoa_rise) {
            return ((toa << 7) - (ftoa_rise << 3) + (uftoa_start - uftoa_stop));
        }

        // Corrects latency delay due to DDLL clock distribution. Units are period of 40MHz (25ns)
        uint64_t toa_clkdll_correction(uint64_t spgroup_addr = 0) {
            return spgroup_addr << 2;
        }

        // address including pixel, super pixel and super pixel group values
        uint64_t getAddr(uint64_t packet) { return (packet >> 46) & 0x3ffff; }

        // super pixel group address
        uint64_t getSuperPixelGroup(uint64_t packet) { return (packet >> 51) & 0xf; }

        // super pixel address
        uint64_t getSuperPixel(uint64_t packet) { return (packet >> 49) & 0x3; }

        // pixel address
        uint64_t getPixel(uint64_t packet) { return (packet >> 46) & 0x7; }

        // Time of Arrival (ToA) | units of 25 ns (1/(40 MHz))
        uint16_t getToA(uint64_t packet) { return (packet >> 30) & 0xffff; }

        // fine ToA rising edge | units of ~1.56 ns (1/(640 MHz)
        uint64_t getFToARise(uint64_t packet) { return (packet >> 17) & 0x1f; }

        // fine ToA falling edge | units of ~1.56 ns (1/(640 MHz)
        uint64_t getFToAFall(uint64_t packet) { return (packet >> 12) & 0x1f; }

        // Time over Threshold | units of 25 ns  (1/(40 MHz))
        uint64_t getToT(uint64_t packet) { return (packet >> 1) & 0x7ff; }

        // bit to register whether another hit was coming in while pixel was busy
        uint64_t getPileUp(uint64_t packet) { return packet & 0x1; }

        //========================================
        // uftoa encoding change to actual values (original, TPX4TOOLS by kevin heijhoff)
        // =======================================
        uint64_t uftoaBin[16] = {4, 5, 8, 6, 8, 8, 8, 7, 3, 8, 8, 8, 2, 8, 1, 0};
        uint64_t UftoaStart(uint64_t value) { return uftoaBin[value >> 26 & 0x000F]; }
        uint64_t UftoaStop(uint64_t value) { return uftoaBin[value >> 22 & 0x000F]; }
    };
} // namespace corryvreckan
#endif // Timepix4EVENTLOADER_H
