# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

CORRYVRECKAN_ENABLE_DEFAULT(OFF)
# Define module and return the generated name as MODULE_NAME
CORRYVRECKAN_DUT_MODULE(MODULE_NAME)

# Add source files to library
FIND_PACKAGE(eudaq 2.4 REQUIRED NO_CMAKE_PACKAGE_REGISTRY NO_CMAKE_SYSTEM_PATH)

CORRYVRECKAN_MODULE_SOURCES(${MODULE_NAME} AnalysisItkStripEfficiency.cpp)
TARGET_LINK_LIBRARIES(${MODULE_NAME} ${EUDAQ2_LIBRARY} eudaq::core)

# Provide standard install target
CORRYVRECKAN_MODULE_INSTALL(${MODULE_NAME})