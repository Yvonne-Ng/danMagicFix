#!/bin/sh -v
#PBS -l cput=31:59:59

# Leave combinations of three capital letters as they are, these are replaced during batch submission script i.e. Run_Limits.py

# Direct output
#PBS -k eo
#PBS -e EEE_e.txt
#PBS -o OOO_o.txt

# Go to temporary directory
cd ${TMPDIR}

# Copy inputs
cp -r YYY/* .  

# Setup stuff
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
rcSetup Base,2.3.23 # Careful make match version you're using!

cd Bayesian
#echo "sourcing Setup.sh"
#source Setup.sh

# Run
ZZZ
  
