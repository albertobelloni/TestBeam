#!/bin/bash

###############################################################################
#
# example submission on HEPCMS condor:
#
#    condor_submit condor_fits.jdl
#
###############################################################################

if [ `hostname -d` = "cern.ch" ]; then
    echo "I did not manage to get this to run on LXPLUS, move to HEPCMS"
    exit
fi

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

tiles=('EJ_260' 'EJ_260_2P' 'EJ_200' 'SCSN_81F1' \
    'SCSN_81F2' 'SCSN_81F3' 'SCSN_81F4' 'SCSN_81S')

root -b -q run_fitter.C\(\"energy_tree_${tiles[$1]}\"\)

mkdir -p condorfits_$2

if [ -e $1.stdout ] && [ -e $1.stderr ] && [ -e $1.log ]; then
    mv $1.stdout $1.stderr $1.log condorfits_$2/
fi

exit

