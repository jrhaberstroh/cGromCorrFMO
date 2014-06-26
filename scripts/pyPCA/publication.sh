#!/bin/bash

dir=~/Dropbox/Physics/subgroup/2014-06-13/2ns_data
mkdir $dir

for i in {1..7}; do
    python f2_1TrackModes.py $i -modes1dhist 1 2 3 -1dparabola -plotspectrum -savemode pdf png -outfnamebase $dir/plot -modes2dhist 1 2 -modes2dhist 2 3 -modes2dhist 1 3 -modesCt 1 2 3
done
