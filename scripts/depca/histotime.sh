#!/bin/bash

for t_ns in {0..49}; do
    python f2_1TrackModes.py -offset $(echo $t_ns)00000 -Nframes 100000 -hist1d 1 2 3
done
