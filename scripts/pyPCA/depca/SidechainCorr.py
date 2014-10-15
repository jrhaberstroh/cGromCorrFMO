import numpy as np
import h5py
import matplotlib.pyplot as plt


def AvgAndCorrelateSidechains(E_t_ia, fEnd = 3, fStart = 0, fStride=1):
        assert(fStart < fEnd)
	numFrames = min(fEnd-fStart, E_t_ia.shape[0]-fStart) / fStride
        fEnd = fStart + (numFrames * fStride)
   
        num_chromo = E_t_ia.shape[1]
        num_vars   = E_t_ia.shape[2]
        Corr_i_ab = np.zeros( (num_chromo, num_vars, num_vars))
        AvgEia  = np.zeros( (num_chromo, num_vars) )
        
        
        for i in xrange(num_chromo):
            print "Computing covariance on chromophore number {}".format(i+1)
            print "\tArray size: ({},{})".format(E_t_ia.shape[2], E_t_ia.shape[0])
            max_floats = 4E9 / 8
            max_times = max_floats / E_t_ia.shape[2]
            dset_times = (fEnd-fStart) / fStride
            chunks = int(np.ceil(dset_times / max_times))
            chunk_size = dset_times/chunks
            print("\tAssuming 4GB RAM usage...\n" +  
                  "\tDesired number of chunks: {}, Chunk size: {} [{}MB]\n".format(chunks, chunk_size, chunk_size*8*E_t_ia.shape[2]/1E8) + 
                  "\tMissing datapoints due to chunk truncation: {}".format((fEnd - fStart) - (chunk_size * chunks * fStride)))
            if (chunks > 1):
                raise NotImplementedError("No chunk feature yet implemented")
            print "Loading data for chromophore {}...".format(i)
            RAM_Datasubset = E_t_ia[fStart:fEnd:fStride,i,:]
            print "Computing covariance for chromophore {}...".format(i)
            Corr_i_ab[i,:,:] = np.cov(RAM_Datasubset, rowvar=0 )
            print "Computing mean for chromophore {}...".format(i)
            AvgEia[i,:]  = RAM_Datasubset.sum(axis=0)
            AvgEia[i,:]  /= numFrames

	return Corr_i_ab, AvgEia


def TimeCorr(f1, f2):
    assert f1.size == f2.size
    if any(f1 != 0.0) or any(f2 != 0.0):
        corrarr = np.correlate(f1, f2, 'full')[f1.size-1:]
        numavg  = corrarr.size - np.arange(corrarr.size)
        corrarr /= numavg
        return corrarr

def hdfSpectralCorr(hdf_file, timetag):
    site = 0
    sc1 = 3
    sc2 = 8
    offset = 0
    nframes = 50000
    start = offset
    stop  = offset + nframes
    with h5py.File(hdf_file, 'r') as f:
        E_tij = f[timetag]
        dt = E_tij.attrs['dt']
        print E_tij.shape
        f1_t = E_tij[start:stop,site,sc1] - np.mean(E_tij[:,site,sc1])
        f2_t = E_tij[start:stop,site,sc2] - np.mean(E_tij[:,site,sc2])
        mycorr = TimeCorr(f1_t, f2_t)
        mycorr = mycorr[0:len(mycorr)/10]
# SPECTRAL DENSITY CODE
#        Jw = np.absolute(np.fft.rfft(mycorr))
#        w_cm  = np.fft.rfftfreq(len(mycorr),dt) * 5.308 * 2. * np.pi
#        print w_cm.shape
#        print Jw.shape
#        #plt.plot(w_cm, Jw * w_cm)
#        N = 16
#        smooth = pd.rolling_mean(Jw * w_cm, N)
#        smooth = np.roll(smooth, -N/2)
#        plt.plot(w_cm, smooth, linewidth=2)
#        plt.show()
    

def HDFStoreCt_v1(hdf_file, timetag, cttag, site_a, site_b, chromo=1, offset = 0, nframes = 200000):
    start = offset
    stop  = offset + nframes
    chromo -= 1
    with h5py.File(hdf_file, 'r') as f:
        E_tij = f[timetag]
        dt = E_tij.attrs['dt']
        print E_tij.shape
        f1_t = E_tij[start:stop,chromo,site_a] - np.mean(E_tij[:,chromo,site_a])
        f2_t = E_tij[start:stop,chromo,site_b] - np.mean(E_tij[:,chromo,site_b])
    mycorr = TimeCorr(f1_t, f2_t)
    mycorr = mycorr[0:len(mycorr)/10]
    return mycorr
    #plt.plot(mycorr)
    #plt.show()
    
    #with h5py.File(hdf_file, 'r') as f:
    #    C_abt = f[cttag]


def HDFStoreCt_v2(hdf_file, timetag, cttag, site_a, site_b, chromo=1, offset = 0, nframes = 200000, corr_len=20000):
    start = offset
    stop  = offset + nframes
    chromo -= 1
    mycorr = np.zeros((corr_len))

    f1_t = np.zeros((stop-start))
    f2_t = np.zeros((stop-start))
    #lock.acquire()
    try:
        with h5py.File(hdf_file, 'r') as f:
            E_tij = f[timetag]
            dt = E_tij.attrs['dt']
            print E_tij.shape
            if (start >= E_tij.shape[0]):
                raise RuntimeError("offset > len(E_tij) [{} > {}]; Use a smaller offset.".format(start,E_tij.shape[0]))
            if (stop > E_tij.shape[0]):
                raise RuntimeError("offset + nframes > len(E_tij) [{} > {}]; Use a smaller offset or fewer frames.".format(stop,E_tij.shape[0]))
            f1_t[:] = E_tij[start:stop,chromo,site_a] - np.mean(E_tij[:,chromo,site_a])
            f2_t[:] = E_tij[start:stop,chromo,site_b] - np.mean(E_tij[:,chromo,site_b])
    except RuntimeError as re:
        #lock.release()
        raise re

    #lock.release()
    #mycorr = TimeCorr(f1_t, f2_t)
    for tc in xrange(corr_len):
        mycorr[tc] = np.inner(f1_t[0:nframes-tc],f2_t[tc:nframes]) / (nframes - tc)
    return mycorr


def HDFStoreCt_v3(hdf_file, timetag, cttag, site_a, site_b, chromo=1, offset = 0, nframes = 200000, corr_len=20000):
    start = offset
    stop  = offset + nframes
    chromo -= 1

    f1_t = np.zeros((stop-start))
    f2_t = np.zeros((stop-start))
    with h5py.File(hdf_file, 'r') as f:
        E_tij = f[timetag]
        dt = E_tij.attrs['dt']
        print E_tij.shape
        if (start >= E_tij.shape[0]):
            raise RuntimeError("offset > len(E_tij) [{} > {}]; Use a smaller offset.".format(start,E_tij.shape[0]))
        if (stop > E_tij.shape[0]):
            raise RuntimeError("offset + nframes > len(E_tij) [{} > {}]; Use a smaller offset or fewer frames.".format(stop,E_tij.shape[0]))
        f1_t[:] = E_tij[start:stop,chromo,site_a] - np.mean(E_tij[:,chromo,site_a])
        f2_t[:] = E_tij[start:stop,chromo,site_b] - np.mean(E_tij[:,chromo,site_b])

    numavg  = nframes - np.arange(corr_len)
    mycorr = np.real(np.fft.ifft( np.fft.fft(f2_t) * np.conj(np.fft.fft(f1_t)) )[0:corr_len])/numavg

    return mycorr

