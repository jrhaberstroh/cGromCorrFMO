#!/bin/bash

style="data"
data="site1_SHOinput.dat"
hdfout=gendata/data/FMO_site1.hdf5
if [ $style == "brownian" ]; then
    echo 'autogendata.sh: Computing brownian dynamics'
    amppeak=100
    pow=2
    hdfout=gendata/data/gen-brn-A$amppeak-pow$pow.hdf5
    python gendata/genSHO.py -dynamics brownian -pow $pow -amppeak $amppeak -h5out $hdfout
elif [ $style == "data" ]; then
    echo 'autogendata.sh: Computing brownian dynamics using FMO data'
    cat $data
    hdfout=gendata/data/FMO_site1.hdf5
    python gendata/genSHO.py -dynamics brownian -setSHO $(cat $data) -h5out $hdfout -T 20000 -dt .1
elif [ $style == "micro" ]; then
    echo 'autogendata.sh: Computing microcanonical dynamics'
    DOF=9
    amppeak=100
    pow=3
    hdfout=gendata/data/gen-A$amppeak-pow$pow-DOF$DOF.hdf5
    python gendata/genSHO.py -DOF $DOF -pow $pow -amppeak $amppeak -h5out $hdfout
fi

echo 'autogendata.sh: Computing PCA'
cd pyPCA
sed -i "s#^h5file = .*#h5file = ../$hdfout#" f0postProcess.cfg
python f1SidechainCorr.py .01
python f2_1TrackModes.py 1 -plotspectrum -modes1dhist 1 2 3 -1dparabola

