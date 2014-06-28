cGromCorr README
------

Package contains:
    /src, containing cGromCorr executable. This takes an FMO trajectory (.gro) and computes the CDC interactions (Renger & Madjet 2006) to get the dynamics and decomposition of the energy gap in time.
    /scripts/pyPCA, a package of executables to process the output from cGromCorr and perform analysis on the results. As of 4/14, this includes PCA, mode timeseries, histograms, and spatial visualization of modes via Chimera.
    /scripts/gendata, a package of python scripts to generate SHO data in the same format to compare to the outputs of simulation data

Requirements:
    * boost_filesystem (reqd: Boost v1.55)
    * boost_system (reqd: Boost v1.55)
    * c++11 (reqd: gcc v4.8.0)


dE2csv utilizes:
    * A reduced "forcefield" (charges and masses), to be stored in a dictionary indexed by (molName,atomNames). 
        Data is loaded from a .top file with .itp includes.
    * Data from .gro file to locate atomNames (required to have chromophores called BCL), loaded frame-by-frame.
    * A .itp file specifying the chromophore "CDC" ground/excited state charges (required to be called BCX).


Basic run instructions:
1. Make using src/basic_test.sh for basic compilation check
2. Locate and connect FMO topology file (BCL atomname required)
3. Select a number of frames to read that is <= number of frames in the file, for passing into basic_run.sh 
(WARNING: WILL SEG FAULT IF CONDITION IT NOT MET.)
4. Make and run using src/basic_run.sh 
