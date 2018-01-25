#!bin/bash

rm -rf build/
rm -rf install/
mkdir build/
cd build/
cmake ..
cp ../make.sh .
source make.sh

