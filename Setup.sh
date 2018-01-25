#!bin/bash

setupATLAS
#asetup AtlasDerivation,21.0.19.1,here
#asetup AnalysisBaseExternals-21.2.0,here
lsetup "asetup AnalysisBase,21.2.3,here"
#lsetup "asetup AnalysisBase,21.1,here"
lsetup root
export PYTHONPATH=./install/InstallArea/x86_64-slc6-gcc62-opt/src/Bayesian/plotting/PythonModules:${PYTHONPATH}

#setup numpy
#echo
#echo "Running pyanalyssi setuppppp"
#lsetup /cvmfs/sft.cern.ch/lcg/releases/LCG_88/pyanalysis/2.0/x86_64-slc6-gcc62-opt/
#export CPP_INCLUDE_PATH=/lustre/SCRATCH/atlas/ywng/r21/r21Rebuild/r21StatisticalAnalysis/yaml-cpp/install/usr/local/include:$CPP_INCLUDE_PATH
#export LD_LIBRARY_PATH=/lustre/SCRATCH/atlas/ywng/r21/r21Rebuild/r21StatisticalAnalysis/yaml-cpp/install/usr/local/lib:$LD_LIBRARY_PATH
#localSetupSFT /cvmfs/sft.cern.ch/lcg/releases/LCG_88/pyanalysis/2.0/x86_64-slc6-gcc62-opt/
#localSetupSFT pyanalys#is 
#/cvmfs/atlas.cern.ch/r#epo/sw/software/21.2/sw/lcg/releases/pyanalysis/2.0-32412
cd install/InstallArea/x86_64-slc6-gcc62-opt
source setup.sh

cd ../../../

