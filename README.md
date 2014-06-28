cGromCorr README
=======

##Package contents

* **/src**, containing cGromCorr executable. This takes an FMO trajectory (.gro) and computes the CDC interactions (Renger & Madjet 2006) to get the dynamics and decomposition of the energy gap in time.
* **/scripts/pyPCA**, a package of executables to process the output from cGromCorr and perform analysis on the results. As of 4/14, this includes PCA, mode timeseries, histograms, and spatial visualization of modes via Chimera.
* **/scripts/gendata**, a package of python scripts to generate SHO data in the same format to compare to the outputs of simulation data

##Requirements
* boost_filesystem (req == Boost v1.55)
* boost_system (req == Boost v1.55)
* c++11  (req >= gcc v4.8.0) 


##Basic run instructions
1. Edit GNUmakefile to reflect the location of boost_filesystem and boost_system.
1. Make using src/basic_test.sh for unit tests
1. Edit src/basic_run
  1. Select a trajcectory, topology file, and excited-state forcefield (see API below for requirements)
  1. Select a number of frames to read. (See API for warnings and details)
1. Make and run using src/basic_run.sh 


###API
./traj2dEcsv traj_in topology excited_itp csv_out num_frames
* **traj_in**: Trajectory to process, in gromacs text format .gro. Used to generate atomNames. Molecules with the name BCL are designated as chromophores. File is loaded frame-by-frame, so there are no file size restrictions beyond system restrictions. 
* **topology**: Gromacs format forcefield for the non-chromophore sites (i.e. the environment). Only collects data from charges and masses, stored in a dictionary indexed by (molName,atomName). #include statements will are treated.
* **excited_itp**: A .itp file specifying the chromophore "CDC" ground/excited state charges. The data in the itp file must specify ground state charges and excited state charges as A/B states (in the Gromacs format), and the molecule must be named BCX.
* **csv_out**: Location to save csv data out. Data is saved to csv_out.time and csv_out.mtx for the timeseries and the correlation matrix respectively.
* **num_frames**: Number of frames to analyze for the timeseries. WARNING: if this number is greater than the number of frames in traj_in, program will seg-fault at the EOF. However, this failure should not cause the timeseries to suffer.
