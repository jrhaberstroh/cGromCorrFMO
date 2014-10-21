import h5py
import ConfigParser
import sys
import iopro

def dEhdf5_init(hdf_file, hdf_dsname, open_flag, Ncoarse=0, ntimes=0, dt_ps=None, chunktime=1000):
    with h5py.File(hdf_file,open_flag) as h5_out:
        if not hdf_dsname:
            return
        # Check if the dataset exists
        h5keys = h5_out.items()
        goodkeys = [key[0] == hdf_dsname for key in h5keys]
        if any(goodkeys):
            ds = h5_out[hdf_dsname]
        else:
            ds = h5_out.create_dataset(hdf_dsname, shape=(ntimes,Ncoarse), chunks=(chunktime,Ncoarse), maxshape=(None, Ncoarse), dtype='f32')

        if dt_ps:
            ds.attrs['dt_unit'] = "picoseconds"
            ds.attrs['dt'] = dt_ps

# Python recipe 577096
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


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
            self.Nsites = config.getint('sidechain','Nsites')
            self.csv_files = config.get('sidechain','csvfiles')

            self.sc_file = None
            self.stat_file = None
            self.pca_file= None

    def __enter__(self):
        self.sc_file = None
        self.stat_file = None
        self.pca_file= None
        return self

    def InitSidechain_hdf(self):
        print "Initializing Sidechain data..."
        if not query_yes_no("Are you sure you want to re-write {}?".format(self.sc_h5file), default="no"):
            print "File rewrite skipped."
            return
        Nsites = self.Nsites

        # Load all of the CSV file paths into an array
        with open(self.csv_files) as f:
            csv_files = f.readlines()
        print("Loading data from {} csv files...".format(len(csv_files)))
        
        dset_tags = [self.time_h5tag + "%d" % i for i in xrange(1,Nsites+1)]

        # Load each CSV sequentially into an HDF file that is created on the first step
        first_file=True
        for i,file in enumerate(csv_files):
            print "File #{}: {}".format(i, file.strip())
            x = iopro.loadtxt(file.strip(),delimiter=',')
            assert(format(x[0::Nsites].shape ==  x[(Nsites-1)::Nsites].shape))
            assert(len(x) % Nsites == 0)
            print "\tX total shape: {}".format((x.shape[0]/Nsites, Nsites, len(x[0,:])))
            Ntimes = len(x) / Nsites

            # Create the HDF file if we are looping through for the first time using the sizes from the csv files
            if first_file:
                Ncoarse = len(x[0])
                print "\tCreating new datasets, loaded file with Ncoarse = {} and Ntimes = {}".format(Ncoarse, Ntimes)
                dEhdf5_init(self.sc_h5file, None, 'w')
                for dset_tag in dset_tags:
                    dEhdf5_init(self.sc_h5file, dset_tag, 'a', Ncoarse=Ncoarse, dt_ps=.005)
                first_file=False

            h5_out =  h5py.File(self.sc_h5file,'a')
            try:
                dsets  = [h5_out[dset_tag] for dset_tag in dset_tags]
                for i,ds in enumerate(dsets):
                    oldlen = ds.shape[0]
                    newlen = oldlen + Ntimes
                    ds.resize(newlen, axis=0)
                    ds[oldlen:newlen,:] = x[i::Nsites,:]
                h5_out.close()
            except:
                print "Write failed..."
                h5_out.close()
                print("Error: {}".format(sys.exc_info()[0]))
                raise
    def ExamineSidechain_hdf(self):
        should_close = False
        if not self.sc_file:
            self.sc_file = h5py.File(self.sc_h5file)
            should_close = True
        print self.sc_file.keys()
        if should_close:
            self.sc_file.close()
    def GetSidechain_hdf(self, i):
        if not self.sc_file:
            self.sc_file  = h5py.File(self.sc_h5file)
        return self.sc_file[self.time_h5tag + "%d" % i]
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

