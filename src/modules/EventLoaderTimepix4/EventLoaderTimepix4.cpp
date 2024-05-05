/**
 * @file
 * @brief Implementation of module EventLoaderTimepix4
 * For now using code from the Kepler TPX4 data decoder
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EventLoaderTimepix4.h"

#include <bitset>
#include <cmath>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace corryvreckan;
using namespace std;

EventLoaderTimepix4::EventLoaderTimepix4(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), m_detector(detector), m_currentEvent(0), m_prevTime(0), m_shutterOpen(false) {

    config_.setDefault<size_t>("buffer_depth", 10000);
    m_buffer_depth = config_.get<size_t>("buffer_depth");


    // Take input directory from global parameters
    m_inputDirectory = config_.getPath("input_directory");

}

void EventLoaderTimepix4::initialize() {

    if(m_buffer_depth < 1) {
        throw InvalidValueError(config_, "buffer_depth", "Buffer depth must be larger than 0.");
    } else {
        LOG(INFO) << "Using buffer_depth = " << m_buffer_depth;
    }

    // File structure is RunX/ChipID/files.dat

    // Open the root directory
    DIR* directory = opendir(m_inputDirectory.c_str());
    if(directory == nullptr) {
        throw ModuleError("Directory " + m_inputDirectory + " does not exist");
    } else {
        LOG(TRACE) << "Found directory " << m_inputDirectory;
    }
    dirent* entry;
    dirent* file;

    // Buffer for file names:
    std::vector<std::string> detector_files;
    // Read the entries in the folder
    while((entry = readdir(directory))) {

        // Ignore UNIX functional directories:
        if(std::string(entry->d_name).at(0) == '.') {
            continue;
        }

        // If these are folders then the name is the chip ID
        // For some dt_type is empty for the data entries so currently identifying files by name
        if(entry->d_type == DT_DIR || (entry->d_type == DT_UNKNOWN && std::string(entry->d_name).at(0) == 'T')) {

            // Only read files from correct directory: Note: currently doesn't work? maybe some weird metadata that determines this that isn't set in the tpx4 data file unlike the tpx3 data?
            if(entry->d_name != m_detector->getName()) {
                continue;
            }
            LOG(DEBUG) << "Found directory for detector " << entry->d_name;

            // Open the folder for this device
            string dataDirName = m_inputDirectory + "/" + entry->d_name;
            DIR* dataDir = opendir(dataDirName.c_str());

            // Get all of the files for this chip
            while((file = readdir(dataDir))) {
                string filename = dataDirName + "/" + file->d_name;

                // Check if file has extension .dat
                if(string(file->d_name).find(".dat") != string::npos) {
                    LOG(INFO) << "Enqueuing data file for " << entry->d_name << ": " << filename;
                    detector_files.push_back(filename);
                }
            }
        }
    }

    // Check that we have files for this detector and sort them correctly:
    if(detector_files.empty()) {
        throw ModuleError("No data file found for detector " + m_detector->getName() + " in input directory " +
                          m_inputDirectory);
    }

    // Initialise null values for later
    m_syncTime = 0;
    m_clearedHeader = false;
    m_syncTimeTDC = 0;
    m_TDCoverflowCounter = 0;
    m_prevTriggerNumber = 0;
    m_triggerOverflowCounter = 0;
    eof_reached = false;

    // Sort all files by extracting the "serial number" from the file name while ignoring the timestamp:
    std::sort(detector_files.begin(), detector_files.end(), [](std::string a, std::string b) {
        auto get_serial = [](std::string name) {
            const auto pos1 = name.find_last_of('-');
            const auto pos2 = name.find_last_of('.');
            return name.substr(pos1 + 1, pos2 - pos1 - 1);
        };
        return std::stoi(get_serial(a)) < std::stoi(get_serial(b));
    });
    // Open them:
    for(auto& filename : detector_files) {

        auto new_file = std::make_unique<std::ifstream>(filename);
        if(new_file->is_open()) {
            LOG(INFO) << "Opened data file for " << m_detector->getName() << ": " << filename;

            // The header is repeated in every new data file, thus skip it for all.
            uint64_t headerID;

            if(!new_file->read(reinterpret_cast<char*>(&headerID), sizeof headerID)) {
                  throw ModuleError("Cannot read header ID for " + m_detector->getName() + " in file " + filename);
            }
            LOG(WARNING) << "Header ID: \"" << headerID << "\"";
            // comparing 8 byte header with headerID of the file, should be SPIDR4\0\0 which is flipped into 4RDIPS which is 57527937618003
            if(headerID != 57527937618003) {
                throw ModuleError("Incorrect header ID for " + m_detector->getName() + " in file " + filename + ": " +
                                  std::to_string(headerID));
            }

            // Store the file in the data vector:
            m_files.push_back(std::move(new_file));
        } else {
            throw ModuleError("Could not open data file " + filename);
        }
    }

    // Set the file iterator to the first file for every detector:
    m_file_iterator = m_files.begin();

    // Hit time (in s)
    hHitTime = new TH1F("hitTime", "hitTime", 1000, -0.5, 999.5);

    // ToA
    hRawToA = new TH1F("RawToA", "RawToA", 1 << 16, -0.5, (1 << 16)-0.5);
    hRawExtendedToA = new TH1F("RawExtendedToA", "RawExtendedToA", 1000, 0, 1E9);
    hRawFullToA = new TH1F("RawFullToA", "RawFullToA", 1000, 0, 1E9);
    // ToT
    hRawToT = new TH1F("rawToT", "rawToT", 1000, -0.5, 999.5);
    hRawFullToT = new TH1F("rawFullToT", "rawFullToT", 1000, -0.5, 999.5);
    hToT = new TH1F("ToT", "ToT", 1000, -0.5E-9, 999.5E-9);
    // Make debugging plots
    std::string title = m_detector->getName() + " Hit map;x [px];y [px];# entries";
    hHitMap = new TH2F("hitMap",
                       title.c_str(),
                       m_detector->nPixels().X(),
                       -0.5,
                       m_detector->nPixels().X() - 0.5,
                       m_detector->nPixels().Y(),
                       -0.5,
                       m_detector->nPixels().Y() - 0.5);
}

StatusCode EventLoaderTimepix4::run(const std::shared_ptr<Clipboard>& clipboard) {

    // This will loop through each timepix4 registered, and load data from each of them. This can
    // be done in one of two ways: by taking all data in the time interval (t,t+delta), or by
    // loading a fixed number of pixels (ie. 2000 at a time)

    // Check if event frame is defined:
    auto event = clipboard->getEvent();

    LOG(TRACE) << "== New event";

    // If all files for this detector have been read, end the run:
    if(m_file_iterator == m_files.end()) {
        return StatusCode::Failure;
    }

    // Make a new container for the data
    PixelVector deviceData;
    SpidrSignalVector spidrData;

    // Load the next chunk of data
    bool data = loadData(clipboard, deviceData, spidrData);

    // If data was loaded then put it on the clipboard
    if(data) {
        LOG(DEBUG) << "Loaded " << deviceData.size() << " pixels for device " << m_detector->getName();
        clipboard->putData(deviceData, m_detector->getName());
    }

    if(!spidrData.empty()) {
        clipboard->putData(spidrData, m_detector->getName());
    }

    // Otherwise tell event loop to keep running
    LOG_PROGRESS(DEBUG, "tpx3_loader") << "Current time: " << Units::display(event->start(), {"s", "ms", "us", "ns"});

    return StatusCode::Success;
}


bool EventLoaderTimepix4::decodeNextWord() {
    std::string detectorID = m_detector->getName();



    // Clearing databuffer off data from previous word.
    m_dataBuffer.clear();

    LOG(DEBUG) << "Starting word decoding";
    // Check if current file is at its end and move to other one
    if((*m_file_iterator)->eof()) {
        if (!m_file_index){
            m_file_iterator++;
        }
        else{
            m_file_iterator--;
        }
        // If that file also has reached eof then stop here
        if ((*m_file_iterator)->eof()){
            LOG(INFO) << "EOF for all files of " << detectorID;
            eof_reached = true;
            return false;
        }
        LOG(INFO) << "Starting to read next file for " << detectorID << ": " << (*m_file_iterator).get();
    }

    // Now read the data packets.
    uint64_t header = 0;

    // If we can't read data anymore, jump to begin of loop:
    if(!(*m_file_iterator)->read(reinterpret_cast<char*>(&header), sizeof header)) {
        LOG(INFO) << "No more data in current file for " << detectorID << ": " << (*m_file_iterator).get();
        return true;
    }

    LOG(TRACE) << "0x" << hex << header << dec << " - " << header;
    auto [groupID, encoding, contentID, streamID, contentSize] = decode_header(header);
    LOG(DEBUG) << "Group ID " << groupID;
    LOG(DEBUG) << "Content Encoding " << encoding;
    LOG(DEBUG) << "Content ID " << contentID;
    LOG(DEBUG) << "Stream ID " << streamID;
    LOG(DEBUG) << "Content size " << contentSize;
    // Group id 0x7 is undefined non SPIDR user data added via RC
    contentSize = contentSize << 3; // multiplying by 8
    if (m_dataBuffer.size()*8 < contentSize)
        m_dataBuffer.resize(contentSize/8);

    if (groupID == 0x7){
        LOG(TRACE) << "Found user defined data";
        if(!(*m_file_iterator)->read(reinterpret_cast<char*>(m_dataBuffer.data()), contentSize)){
            LOG(INFO) << "No more data in current file for " << detectorID << ": " << (*m_file_iterator).get();
            return true;
        }
        LOG(TRACE) << "User information, skipping!";
    }
    else if (groupID == 0x0){
        LOG(DEBUG) << "Found timepix4 data";
        if (encoding == 0b00){
            if(!(*m_file_iterator)->read(reinterpret_cast<char*>(m_dataBuffer.data()), contentSize)){
                LOG(INFO) << "No more data in current file for " << detectorID << ": " << (*m_file_iterator).get();
                return true;
            }
            LOG(DEBUG) << "Found " << m_dataBuffer.size() << " data packets." << std::endl;
            for (ulong i = 0; i < m_dataBuffer.size(); i++){
                m_dataPacket = m_dataBuffer[i];
                if (decodePacket(m_dataPacket, 16, true)){
                    LOG(TRACE) << "Found pixel data!";
                    bool digCompare = false;
                    for (const auto &digColRow: m_digColRow){
                        digCompare = compareTupleEq(digColRow, m_colrow);
                        if (digCompare) break;
                    }
                    if (!digCompare){
                        uint32_t col = std::get<0>(m_colrow);
                        uint32_t row = std::get<1>(m_colrow);
                        double correctedTime = static_cast<double>(m_fullToa) * 1/(8*640E6); // time in s
                        double correctedToT = static_cast<double>(m_fullTot) * 1/(8*640E6); // tot in s
                        auto pixel = std::make_shared<Pixel>(detectorID, col, row, static_cast<int>(m_fullTot), correctedToT, m_heartbeat);
                        pixel->setCharge(correctedToT);
                        sorted_pixels_.push(pixel);
                        hRawToT->Fill(m_tot);
                        hRawFullToT->Fill(m_fullTot);
                        hToT->Fill(correctedToT);
                        hHitMap->Fill(col, row);
                        hRawToA->Fill(m_toa);
                        hRawExtendedToA->Fill(m_ext_toa);
                        hRawFullToA->Fill(m_fullToa);
                        hHitTime->Fill(correctedTime);

                    }
                }
                else{
                    LOG(TRACE) << "Found heartbeat data!";
                }
            }
        }
        else {
            LOG(TRACE) << "Pixel encoding wrong, skipping!";        
        }
    }
    else{
        if(!(*m_file_iterator)->read(reinterpret_cast<char*>(m_dataBuffer.data()), contentSize)){
            LOG(INFO) << "No more data in current file for " << detectorID << ": " << (*m_file_iterator).get();
            return true;
        }
        LOG(WARNING) << "Other type of data, ignored for now";
    }

    // noting down original position within file stream before switching

    m_stream_pos[m_file_index] = (*m_file_iterator)->tellg();
    LOG(DEBUG) << "Current stream index position " << m_stream_pos[m_file_index];
    LOG(DEBUG) << "Finished reading event from file " << m_file_index;
    if (m_file_index == 0){
        m_file_index = 1;
        m_file_iterator++;
    }
    else if (m_file_index == 1){
        m_file_index = 0;
        m_file_iterator--;
    }
    else {
        LOG(WARNING) << "File index in non sensical position";
    }
    return true;
}

bool EventLoaderTimepix4::decodePacket(uint64_t dataPacket, uint64_t ratio_VCO_CKDLL=16, bool gray=true){
    uint64_t top = (dataPacket >> 63) & 0x1;
    uint64_t header = (dataPacket >> 55) & 0xFF;
//    if (!top) LOG(WARNING) << "Top? " << top;
    if (header > 0xDF) { // heartbeat data
        m_heartbeat = dataPacket & 0x7FFFFFFFFFFFFFull;
        LOG(TRACE) << "Heartbeat data: " << m_heartbeat;
        return false;
    }
    else { // pixel data
        m_addr = (dataPacket>>46) & 0x3ffff; // address including pixel, super pixel and super pixel group values
        m_sPGroup = (dataPacket >> 51) & 0xf; // super pixel group address
        m_sPixel = (dataPacket >> 49) & 0x3; // super pixel address
        m_pixel = (dataPacket >> 46) & 0x7; // pixel address

        m_toa = (dataPacket >> 30) & 0xffff; // Time of Arrival (ToA) | units of 25 ns (1/(40 MHz))
        if (gray) m_toa = GrayToBin(m_toa);
        m_ftoa_rise = (dataPacket >> 17) & 0x1f; // fine ToA rising edge | units of ~1.56 ns (1/(640 MHz)
        m_ftoa_fall = (dataPacket >> 12) & 0x1f; // fine ToA falling edge | units of ~1.56 ns (1/(640 MHz)
        m_tot = (dataPacket>>1) & 0x7ff; // Time over Threshold | units of 25 ns  (1/(40 MHz))
        m_pileup = dataPacket & 0x1; // bit to register whether another hit was coming in while pixel was busy

        m_uftoa_start = UftoaStart(dataPacket); // ultra fine ToA start encoding | units of ~195 ps (1/(640*8 MHz))
        m_uftoa_stop = UftoaStop(dataPacket); // ultra fine ToA stop encoding | units of ~195 ps (1/(640*8 MHz))


        m_ext_toa = extendToa(m_toa, m_heartbeat, m_tot);
        m_fullTot=fullTot(m_ftoa_rise, m_ftoa_fall, m_uftoa_start, m_uftoa_stop, m_tot); // full corrected ToT | units of ~195 ps (1/(640*8 MHz))
        m_fullToa=fullToa(m_ext_toa ,m_uftoa_start, m_uftoa_stop, m_ftoa_rise) - toa_clkdll_correction(m_sPGroup); // rfull corrected ToA | units of ~195 ps (1/(640*8 MHz))
        m_colrow = decodeColRow(m_pixel, m_sPixel, m_sPGroup, header, top); // decodes the row and col value from the address

        LOG(TRACE) << " ";
        LOG(TRACE) << "Col " << std::get<0>(m_colrow);
        LOG(TRACE) << "Row " << std::get<1>(m_colrow);
        LOG(TRACE) << "tot " << m_tot;
        LOG(TRACE) << "ftoa_fall " << m_ftoa_fall;
        LOG(TRACE) << "ftoa_rise " << m_ftoa_rise;
        LOG(TRACE) << "uftoa_start " << m_uftoa_start;
        LOG(TRACE) << "uftoa_stop " << m_uftoa_stop;
        LOG(TRACE) << "toa " << m_toa;
        LOG(TRACE) << "pixel " << m_pixel;
        LOG(TRACE) << "super pixel " << m_sPixel;
        LOG(TRACE) << "fullTot " << m_fullTot;
        LOG(TRACE) << "fullToa " << m_fullToa;
        LOG(TRACE) << "Super Pixel group " << m_sPGroup;
    }
    return true;

}



void EventLoaderTimepix4::fillBuffer() {
    // read data from file and fill timesorted buffer
    while(sorted_pixels_.size() < m_buffer_depth && !eof_reached) {
        // decodeNextWord returns false when EOF is reached and true otherwise
        if(!decodeNextWord()) {
            LOG(TRACE) << "decodeNextWord returns false: reached EOF.";
            break;
        }
    }
}

// Function to load data for a given device, into the relevant container
bool EventLoaderTimepix4::loadData(const std::shared_ptr<Clipboard>& clipboard,
                                   PixelVector& devicedata,
                                   SpidrSignalVector& spidrData) {

    std::string detectorID = m_detector->getName();
    auto event = clipboard->getEvent();

    LOG(DEBUG) << "Loading data for device " << detectorID;
    fillBuffer();

    // Now we have data buffered into the temporary storage. We will sort this by time, and then load
    // the data from one event onto it.

    while(!sorted_pixels_.empty()) {
        auto pixel = sorted_pixels_.top();

        auto position = event->getTimestampPosition(pixel->timestamp());

        if(position == Event::Position::AFTER) {
            LOG(DEBUG) << "Stopping processing event, pixel is after "
                          "event window ("
                       << Units::display(pixel->timestamp(), {"s", "us", "ns"}) << " > "
                       << Units::display(event->end(), {"s", "us", "ns"}) << ")";
            break;
        } else if(position == Event::Position::BEFORE) {
            LOG(TRACE) << "Skipping pixel, is before event window (" << Units::display(pixel->timestamp(), {"s", "us", "ns"})
                       << " < " << Units::display(event->start(), {"s", "us", "ns"}) << ")";
            sorted_pixels_.pop();
        } else {
            devicedata.push_back(pixel);
            sorted_pixels_.pop();
        }

        // Refill buffer
        fillBuffer();
    }

    while(!sorted_signals_.empty()) {
        auto signal = sorted_signals_.top();

        auto position = event->getTimestampPosition(signal->timestamp());

        if(position == Event::Position::AFTER) {
            LOG(DEBUG) << "Stopping processing event, signal is after "
                          "event window ("
                       << Units::display(signal->timestamp(), {"s", "us", "ns"}) << " > "
                       << Units::display(event->end(), {"s", "us", "ns"}) << ")";
            break;
        } else if(position == Event::Position::BEFORE) {
            LOG(TRACE) << "Skipping signal, is before event window ("
                       << Units::display(signal->timestamp(), {"s", "us", "ns"}) << " < "
                       << Units::display(event->start(), {"s", "us", "ns"}) << ")";
            sorted_signals_.pop();
        } else {
            spidrData.push_back(signal);
            sorted_signals_.pop();
        }
    }

    // If no data was loaded, return false
    if(devicedata.empty()) {
        return false;
    }

    // Count events:
    m_currentEvent++;
    return true;
}
