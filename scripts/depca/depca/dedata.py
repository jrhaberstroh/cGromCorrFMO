import h5py
import ConfigParser

class dEData():
    def __init__(self, config = './f0postProcess.cfg'):
        with open(config) as fp:
            config = ConfigParser.RawConfigParser()
            config.readfp(fp)
            self.sc_h5file = config.get('sidechain','h5file')
            self.time_h5tag = config.get('sidechain','time_h5tag')
            self.h5stats= config.get('sidechain','h5stats')
            self.h5eavtag = config.get('sidechain','h5eavtag')
            self.h5corrtag = config.get('sidechain','h5crtag')
            self.pca_h5file = config.get('sidechain','pcafile')

    def __enter__(self):
        self.sc_file = None
        self.stat_file = None
        self.pca_file= None
        return self

    def GetSidechain_hdf(self):
        if not self.sc_file:
            self.sc_file  = h5py.File(self.sc_h5file)
        return self.sc_file[self.time_h5tag]
    def CloseSidechain_hdf(self):
        if self.sc_file:
            self.sc_file.close()
            self.sc_file = None

    def GetStats_hdf(self):
        if not self.stat_file:
            self.stat_file = h5py.File(self.h5stats)
        return self.stat_file[self.h5eavtag], self.stat_file[self.h5corrtag]
    def CloseStats_hdf(self):
        if self.stat_file:
            self.stat_file.close()
            self.stat_file = None

    def GetModes_hdf(self):
        if not self.pca_file:
            self.pca_file = h5py.File(self.pca_h5file)
        return self.pca_file[self.time_h5tag]
    def CloseModes_hdf(self):
        if self.pca_file:
            self.pca_file.close()
            self.pca_file = None

    def __exit__(self, type, value, traceback):
        self.CloseSidechain_hdf()
        self.CloseStats_hdf()
        self.CloseModes_hdf()

