#!/bin/bash

dir=~/Dropbox/Physics/subgroup/2014-06-28/50ns-data
mkdir $dir

for i in {1..7}; do
    python f2_1TrackModes.py -site $i -hist1d 1 2 3 -1dparabola -plotspectrum -savemode pdf png -outfnamebase $dir/plot -hist2d 1 2 -hist2d 2 3 -hist2d 1 3 -modesCt 1 2 3
done

