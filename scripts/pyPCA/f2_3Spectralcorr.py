import ConfigParser
import h5py
import numpy as np
import matplotlib.pyplot as plt

def TimeCorr(f1, f2):
    assert f1.size == f2.size
    if any(f1 != 0.0) or any(f2 != 0.0):
        corrarr = np.correlate(f1, f2, 'full')[f1.size-1:]
        numavg  = corrarr.size - np.arange(corrarr.size)
        corrarr /= numavg
        corrarr /= corrarr[0]
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
        mycorr = mycorr[0:len(mycorr)/2]
        Jw = np.absolute(np.fft.rfft(mycorr))
        w  = np.fft.rfftfreq(len(mycorr),dt)
        print w.shape
        print Jw.shape
        plt.plot(w, Jw * w)
        plt.show()
    

def main():
    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')
    h5_filename = config.get('sidechain','h5file')
    h5time = config.get('sidechain','time_h5tag')
    h5corr = config.get('sidechain','corr_h5tag')

    hdfSpectralCorr(h5_filename,h5time)

if __name__ == "__main__":
    main()
