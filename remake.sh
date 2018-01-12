#!/bin/bash

if [ -f install/InstallArea/x86_64-slc6-gcc62-opt/setup.sh ]
then 
    cd install/InstallArea/x86_64-slc6-gcc62-opt 
    source setup.sh
    echo "sourcing setup done"
    cd ../../..
else 
    echo "the install setup file doesn't exist!"
fi

