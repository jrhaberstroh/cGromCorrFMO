import argparse
import h5py
import numpy as np
import numpy.linalg as LA
import sys
from numbapro import vectorize
from numbapro import guvectorize

def dEhdf5_init(hdf_file, hdf_dsname, open_flag, Ncoarse, ntimes=0, dt_ps=None, Nsite=7, chunktime=5000):
    with h5py.File(hdf_file,open_flag) as h5_out:
        dsshape=(None,Nsite,Ncoarse)
        chunk = (chunktime, 1, Ncoarse)
        # Check if the dataset exists
        h5keys = h5_out.items()
        goodkeys = [key[0] == hdf_dsname for key in h5keys]
        if any(goodkeys):
            ds = h5_out[hdf_dsname]
        else:
            if chunk:
                ds = h5_out.create_dataset(hdf_dsname, shape=(ntimes,dsshape[1],dsshape[2]), chunks=chunk, maxshape=dsshape, dtype='f32')

        if dt_ps:
            ds.attrs['dt_unit'] = "picoseconds"
            ds.attrs['dt'] = dt_ps


def ApplyPCA_hdf5(E_t_ij, E_avg_ij, modes_inj, hdf_file, hdf_dsname, site, create=True, overwrite=False, Nframes=None, offset=0):
	"""
	DeltaModeTracker: Takes E(t) np.ndarray and an iterable of modes, and plots the deviation of dE(t) [after taking the average for each site]
	                  as contributed from the modes, as compared to the actual dE(t)

	\input 	E_t_ij 		- a numpy array of Eij as a function of time
                E_avg_ij        - a numpy array of the mean of E_t_ij
		modes_inj 	- a collection of n normalized modes (in j coordinates) for each site i -- rotation matrix between bases

        \output                 Creates or appends a dataset at hdf_file/hdf_dsname
	"""

	Ntimes  = E_t_ij.shape[0]
	Nsites  = E_t_ij.shape[1]
	Ncoarse = E_t_ij.shape[2]

        # Flip the modes so that they are in order from largest to smallest and zero-based 
	print "Ntimes: ", Ntimes,
	print "Nsites: ", Nsites, 
	print "Ncoarse: ", Ncoarse,

        # Compute t0 and tf from Ntimes and offset
        if Nframes:
            outLen = min(Ntimes-offset,Nframes)
        else:
            outLen = Ntimes
        print outLen
       
        print("Loading output file {}".format(hdf_file))
        open_flag = 'a'
        if overwrite:
            print "\tWARNING: Opening with 'w', file will be overwritten."
            open_flag = 'w'

        with h5py.File(hdf_file, open_flag) as pca_file:
            if create:
                dssize = (outLen, Nsites, Ncoarse)
                pca_tin = pca_file.create_dataset(hdf_dsname, dssize, dtype='f32', chunks=(10000,1,1))
                print("\tCreated Dataset {}, size {}".format(hdf_dsname,dssize))

        t0 = offset
        tf = offset+outLen

	# First, check that the modes are normalized, then weight the modes
        modeweight_inj = np.zeros(modes_inj.shape)
        modeweight_inj[:] = modes_inj[:]
	for m_nj in modeweight_inj:
            for m_j in m_nj:
                if not abs(LA.norm(m_j) - 1) < 1E-5:
                    raise RuntimeError("Bad normalization encountered, {}".format(LA.norm(m_j)-1))
                # Why are my modes squaring??? In the mode matrix???
                m_j *= sum(m_j)
	
    
	#Then, compute mode timeseries:
	print "Computing dE for site {} (O(T*Ncoarse)) using the total mean...".format(site+1)

        GB = float(1E9)
        RAM_GB = .5
        RAM_nfloat = RAM_GB * 1E9 / 8
        RAM_ntimes = RAM_nfloat / Ncoarse
        RAM_nchunk = int( np.ceil((tf - t0) / float(RAM_ntimes)) )
        RAM_time_per_chunk = (tf - t0) / RAM_nchunk
        print "Number of chunks needed: {} of {} GB each".format(RAM_nchunk, RAM_time_per_chunk * Ncoarse * 8 / GB)
        RAM_return_times = (RAM_nchunk*RAM_time_per_chunk)
        RAM_return_tot = (RAM_nchunk*RAM_time_per_chunk) * 8 * (Ncoarse)
        print "Disk space needed for output data: {} GB".format(RAM_return_tot / GB)


        dE_t_j = np.zeros((RAM_time_per_chunk,Ncoarse))
        return_modes_nt = np.zeros((Ncoarse, RAM_return_times))
        dDEresidual_t = np.zeros((RAM_return_times))
        #TODO: Implement chunking to compute large datasets
        for chunk_num in xrange(RAM_nchunk):
            print "Chunk {}:".format(chunk_num+1),;sys.stdout.flush()
            # Each chunk reads [t0_chunk, tf_chunk)
            t0_chunk = (  chunk_num   * RAM_time_per_chunk) + t0
            tf_chunk = ((chunk_num+1) * RAM_time_per_chunk) + t0
            t0_return = t0_chunk - t0
            tf_return = tf_chunk - t0
            print "Loading chunk into RAM...",; sys.stdout.flush()
            RAM_E_t_ij = E_t_ij[t0_chunk:tf_chunk,site,:]

            @vectorize(['float64(float64,float64)'], target='cpu')
            def ParallelSub64(a,b):
                return a - b

            @vectorize(['float32(float32,float32)'], target='cpu')
            def ParallelSub32(a,b):
                return a - b

            # Build dE for chunk
            print "Computing chunk dE...",; sys.stdout.flush()

            try:
                RAM_dE_t_j = ParallelSub32( RAM_E_t_ij, E_avg_ij[site,:])
            except TypeError:
                try:
                    print "32 bit parallel float computation failed, trying 64 bit...",;sys.stdout.flush()
                    RAM_dE_t_j = ParallelSub64( RAM_E_t_ij, E_avg_ij[site,:])
                    print "Success."
                except TypeError as e:
                    print "64 bit failed."
                    raise e
                
            
            print "Rotating chunk...",; sys.stdout.flush()
            RAM_dE_rotated = np.inner(RAM_dE_t_j[:,:], modeweight_inj[site,:,:].T)
        
            print "Writing chunk...",; sys.stdout.flush()
            with h5py.File(hdf_file, 'a') as pca_file:
                pca_tin = pca_file[hdf_dsname]
                pca_tin[t0_chunk:tf_chunk, site, :] = RAM_dE_rotated[:,:]

            print 'Chunk {}  done'.format(chunk_num+1)

	print "Done."
