#!/bin/sh -v

# Leave combinations of three capital letters as they are, these are replaced during batch submission script i.e. Run_Limits.py

# LXBATCH uses a setup already established, so assuming your environment was set up correctly to run, this should work as is.


# Setup stuff
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
cd YYY
source ./rcSetup.sh
#rcSetup Base,2.3.19 # Careful make match version you're using!
#rc find_packages

cd Bayesian

# Run
ZZZ
  
