#!bin/bash

rm -rf build/
mkdir build/
cd build/
cmake ..
cp ../make.sh .
source make.sh

