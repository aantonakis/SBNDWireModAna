#!/bin/sh

SETUP="/cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh"
SBNDCODE_VERSION="v10_06_00"
SBNDCODE_QUAL="e26:prof"
CMAKE_VERSION="v3_27_4"

export SCECORR_PATH="/exp/sbnd/app/users/aantonak/SBND_calib_recom"

source ${SCECORR_PATH}/setup.sh

source ${SETUP}
setup sbndcode ${SBNDCODE_VERSION} -q ${SBNDCODE_QUAL}
setup cmake ${CMAKE_VERSION}
