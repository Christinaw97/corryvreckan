# SPDX-FileCopyrightText: 2017-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

# Option to build all modules
OPTION(BUILD_ALL_MODULES "Build all modules?" OFF)

# reset the saved libraries
SET(CORRYVRECKAN_MODULE_LIBRARIES "" CACHE INTERNAL "Module libraries")

FILE(GLOB subdirs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)
FOREACH(subdir ${subdirs})
    IF(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
        # Internal macros check if the module is actually enabled
        ADD_SUBDIRECTORY(${subdir})
    ENDIF()
ENDFOREACH()
