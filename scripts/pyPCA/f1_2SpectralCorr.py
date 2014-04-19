import ConfigParser
import h5py
import numpy as np

from f1SidechainCorr import h5tstag



def TimeCorr(f1, f2):
    assert f1.size == f2.size
    if any(f1 != 0.0) or any(f2 != 0.0):
        corrarr = np.correlate(f1, f2, 'full')[f1.size-1:]
        numavg  = corrarr.size - np.arange(corrarr.size)
        corrarr /= numavg
        corrarr /= corrarr[0]
        return corrarr


def hdfSpectralCorr(hdf_file, hdf_tag):
    with h5py.File(hdf_file, 'r') as f:
        tsdset = f[hdf_tag+h5tstag]
        print tsdset.shape
        f1 = tsdset[:,0,3] - np.mean(tsdset[:,0,3])
        f2 = tsdset[:,0,8] - np.mean(tsdset[:,0,8])
        mycorr = TimeCorr(f1, f2)
        print mycorr
        print list(mycorr).index(mycorr[1])
    

def main():
    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')
    csv_filename = config.get('sidechain','csv_file')
    h5_filename = config.get('sidechain','h5_file')
    h5_tag = config.get('sidechain','h5_tag')
    t_start = config.getint('sidechain','t_start')
    t_end = config.getint('sidechain','t_end')

    hdfSpectralCorr(h5_filename,h5_tag)

if __name__ == "__main__":
    main()
