#!/bin/bash

#brownian=1
#if [ -z $brownian ]; then
#    echo 'autogendata.sh: Computing microcanonical dynamics'
#    DOF=9
#    amppeak=100
#    pow=3
#    hdfout=gendata/data/gen-A$amppeak-pow$pow-DOF$DOF.hdf5
#    python gendata/genSHO.py -DOF $DOF -pow $pow -amppeak $amppeak -h5out $hdfout
#else
#    echo 'autogendata.sh: Computing brownian dynamics'
#    amppeak=100
#    pow=2
#    hdfout=gendata/data/gen-brn-A$amppeak-pow$pow.hdf5
#    python gendata/genSHO.py -dynamics brownian -pow $pow -amppeak $amppeak -h5out $hdfout
#fi

hdfout="gendata/data/fmo_sho.hdf5"
echo 'autogendata.sh: Computing PCA'
cd pyPCA
sed -i "s#^h5file = .*#h5file = ../$hdfout#" f0postProcess.cfg
python f1SidechainCorr.py .01
python f2_1TrackModes.py 1 -plotspectrum -modes1dhist 1 2 3

