#!/bin/bash

# SPDX-FileCopyrightText: 2022-2024 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

set -e

UPDATE_SCRIPT=$(dirname $0)/update_copyright_years.py
APPENDIX="CERN and the Corryvreckan authors"

python3 ${UPDATE_SCRIPT} -v -a "${APPENDIX}" -x 3rdparty/ -x .reuse/dep5

# update copyright info output corry exe
python3 ${UPDATE_SCRIPT} -v -a "${APPENDIX}" -p "std::cout << \"Copyright (c)" -f src/exec/corry.cpp
