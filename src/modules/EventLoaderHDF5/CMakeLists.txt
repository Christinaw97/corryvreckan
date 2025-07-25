# SPDX-FileCopyrightText: 2023-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

CORRYVRECKAN_ENABLE_DEFAULT(OFF)
CORRYVRECKAN_DETECTOR_TYPE(${MODULE_NAME} "HDF5")

# Define module and return the generated name as MODULE_NAME
CORRYVRECKAN_DETECTOR_MODULE(MODULE_NAME)

# Add source files to library
CORRYVRECKAN_MODULE_SOURCES(${MODULE_NAME} EventLoaderHDF5.cpp)

FIND_PACKAGE(HDF5 COMPONENTS CXX REQUIRED)
TARGET_INCLUDE_DIRECTORIES(${MODULE_NAME} SYSTEM PRIVATE ${HDF5_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${MODULE_NAME} ${HDF5_LIBRARIES})

# Provide standard install target
CORRYVRECKAN_MODULE_INSTALL(${MODULE_NAME})
