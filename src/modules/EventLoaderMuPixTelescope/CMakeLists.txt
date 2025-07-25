# SPDX-FileCopyrightText: 2019-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

CORRYVRECKAN_ENABLE_DEFAULT(OFF)

# Define module and return the generated name as MODULE_NAME
CORRYVRECKAN_GLOBAL_MODULE(MODULE_NAME)

find_package(Boost REQUIRED)

FIND_PACKAGE(MUPIX8DAQ REQUIRED)
find_package(Eigen3 REQUIRED)

IF(NOT MUPIX8DAQ_FOUND)
    MESSAGE(FATAL_ERROR "Could not find MuPix-telescope. You need to do: \n export MUPIX8DAQ=/path/to/repo")
ENDIF()

INCLUDE_DIRECTORIES(SYSTEM ${MUPIX8DAQ_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
message("Setting MUPIX8_DAQ_INCLUDE_DIRS to :${MUPIX8DAQ_INCLUDE_DIRS}")

# Add source files to library
CORRYVRECKAN_MODULE_SOURCES(${MODULE_NAME}
    EventLoaderMuPixTelescope.cpp
    # ADD SOURCE FILES HERE...
)
TARGET_LINK_LIBRARIES(${MODULE_NAME} ${MUPIX8DAQ_LIBRARY})
# Provide standard install target
CORRYVRECKAN_MODULE_INSTALL(${MODULE_NAME})
