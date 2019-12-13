#!/bin/bash

################################################################################
#
# example submission on condor:
#
#    condor_submit condor_fit.jdl
#
# possible values:
#
#    EJ_200, EJ_260, EJ_260_2P, 
#    SCSN_81F1, SCSN_81F2, SCSN_81F3, SCSN_81F4, SCSN_81S
#
# remember to set RUN_DIR and TILE_TO_FIT in the condor_fit.jdl file
#
################################################################################



cd $1
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

root -b -q run_fitter.C\(\"energy_tree_$2\"\)

exit

