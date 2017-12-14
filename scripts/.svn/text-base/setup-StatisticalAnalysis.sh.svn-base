#!/bin/sh

#script to setup the Bayesian Statistical Analysis package

#NOTE to get this script do:
#    svn export svn+ssh://svn.cern.ch/reps/atlasphys-exo/Physics/Exotic/JDM/DiJet/StatisticalAnalysis/Bayesian/trunk/scripts/setup-StatisticalAnalysis.sh

#HOW TO: cd to the directory where you whould like to install the StatisticalAnalysis (Bayesian) code and then source the script; done!

#EXAMPLE:
#    mkdir ~/treasureisland
#    cd ~/treasureisland
#    svn export svn+ssh://svn.cern.ch/reps/atlasphys-exo/Physics/Exotic/JDM/DiJet/StatisticalAnalysis/Bayesian/trunk/scripts/setup-StatisticalAnalysis.sh
#    source setup-StatisticalAnalisis.sh

#NOTE: you may use this script to setup the StatisticalAnalysis code every time you need it, it won't overwrite anything or re-compile unnecessarily


#check where the script is located
#stop if Statistical Analysis package is already available
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
if [[ "$DIR" == *"Bayesian/scripts" ]]; then
    echo -e "\n***WARNING*** running from Bayesian/scripts/: STOP"
    echo -e "***NOTE*** setup using the script located two directories above"
    return 1
fi

#setup RootCore
echo
rcSetup Base,2.3.23
compile=false

#check if Bayesian Statistical Analysis is already present
if [[ ! -d $ROOTCOREBIN/../Bayesian ]]
    then {
        echo -e "\nchecking out Bayesian Statistical Analysis"
	svn co svn+ssh://svn.cern.ch/reps/atlasphys-exo/Physics/Exotic/JDM/DiJet/StatisticalAnalysis/Bayesian/trunk/ Bayesian
	compile=true
    }
else echo -e "\nBayesian Statistical Analysis is already present"
fi

#check if BAT is already present
if [[ ! -d $ROOTCOREBIN/../Asg_BAT ]]
    then {
        echo -e "\nchecking out BAT"
	rc checkout_pkg atlasoff/AsgExternal/Asg_BAT/tags/Asg_BAT-00-09-04-01/ 
	compile=true
    }
else echo -e "\nBAT is already present"
fi

#compile
if [ $compile = true ]; then
    echo -e "\ncompiling"
    rc find_packages
    rc compile    
else
    echo -e "\nno need to compile"
fi

#setup python path
export PYTHONPATH=$ROOTCOREBIN/../Bayesian/plotting/PythonModules/:$PYTHONPATH
#export PYTHONPATH=$PWD/inputs/accDictionary:$PYTHONPATH #is this needed?

#setup numpy
echo
localSetupSFT pyanalysis/1.4_python2.7

echo -e "\ndone"
return 0
