#!/bin/bash

################################################################################
#
# example submission on condor:
#
#    condor_submit condor_fits.jdl
#
#
################################################################################

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

tiles=('EJ_260' 'EJ_260_2P' 'EJ_200' 'SCSN_81F1' \
    'SCSN_81F2' 'SCSN_81F3' 'SCSN_81F4' 'SCSN_81S')

root -b -q run_fitter.C\(\"energy_tree_${tiles[$1]}\"\)

exit

