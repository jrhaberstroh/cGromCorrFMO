import ConfigParser
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import timeit
import multiprocessing as MP

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

## Test code
#    t1 = timeit.Timer('HDFStoreCt_v1("{}","{}", "{}", {}, {})'.format(h5_filename, h5time, h5ct, args.site_a, args.site_b), 'from __main__ import HDFStoreCt_v1').timeit(number=5)
#    print t1
#    t2 = timeit.Timer('HDFStoreCt_v2("{}","{}","{}",{},{})'.format(h5_filename, h5time, h5ct, args.site_a, args.site_b), 'from __main__ import HDFStoreCt_v2').timeit(number=5)
#    print t2
#    t3 = timeit.Timer('HDFStoreCt_v3("{}","{}","{}",{},{})'.format(h5_filename, h5time, h5ct, args.site_a, args.site_b), 'from __main__ import HDFStoreCt_v3').timeit(number=5)
#    print t3
#
#    corr2 = HDFStoreCt_v2(h5_filename, h5time, h5ct, args.site_a, args.site_b)
#    corr3 = HDFStoreCt_v3(h5_filename, h5time, h5ct, args.site_a, args.site_b)
#    corr1 = HDFStoreCt_v1(h5_filename, h5time, h5ct, args.site_a, args.site_b)
#    plt.plot(np.linspace(0,args.lenCt,len(corr2)), corr2)
#    plt.plot(np.linspace(0,args.lenCt,len(corr3)), corr3)
#    plt.plot(np.linspace(0,args.lenCt,len(corr1)), corr1)
#    plt.legend(["Numpy Sums","Fourier Space","Numpy Correlate"])
#    plt.show()

def main():
    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')
    h5_filename = config.get('sidechain','h5file')
    h5time = config.get('sidechain','time_h5tag')
    h5corr = config.get('sidechain','corr_h5tag')
    h5ct = config.get('sidechain','ct_h5tag')

    parser = argparse.ArgumentParser(description = "Module to store and visualize temporal correlations")
    parser.add_argument('-newCt', action="store_true", help="Create new Ct matrix?")
    parser.add_argument('-overwrite', action="store_true", help="Delete old Ct matrix? Must be paired with newCt for confirmation.")
    parser.add_argument('-lenCt', type=float, default=200., help="Length of Ct to store, ps")
    parser.add_argument('-chromo', type=int, default=1, help="Chromophore to compute correlations for")
    args = parser.parse_args()

    num_t = 0
    size_Ct_ab = (0,0,0) # Set value with data from E_tij

    h5ct = "".join([h5ct, "{}".format(args.chromo)])

    with h5py.File(h5_filename, 'r+') as f:
        E_tij = f[h5time]
        num_t = int(np.round(args.lenCt / E_tij.attrs['dt']))
        size_Ct_ab = (num_t, E_tij.shape[2], E_tij.shape[2])

        disk_size = np.product(size_Ct_ab) * 8 / float(10**9)
        print "Note: Estimated maximum size on disk: {0:.1f} GB".format(disk_size)

        if args.newCt and h5ct in f:
            if args.overwrite:
                raise NotImplementedError("No implementation found for overwrite yet!")
            else:
                raise RuntimeError("Warning: Dataset found when newCt was requested, but overwrite flag was not passed. Use -overwrite to delete previous data or remove -newCt flag.")
        elif args.newCt and not h5ct in f:
            f.create_dataset(h5ct, size_Ct_ab)

    print size_Ct_ab
    for a in xrange(size_Ct_ab[1]):
        Ct_b = np.zeros((size_Ct_ab[0],size_Ct_ab[2]))
        for b in xrange(size_Ct_ab[2]):
            print "{},{}".format(a,b),
            Ct_b[:,b] = HDFStoreCt_v3(h5_filename,h5time, h5ct, a, b, chromo=args.chromo)
        with h5py.File(h5_filename, 'r+') as f:
            Ct_ab = f[h5ct]
            print "Data Chunks:",Ct_ab.chunks
            print Ct_ab[:, a, :].shape, Ct_b.shape
            Ct_ab[:,a,:] = Ct_b[:]

if __name__ == "__main__":
    main()
