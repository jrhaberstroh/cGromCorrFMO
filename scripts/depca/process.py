import depca.convert as conv
import depca.sidechain_corr as sc
import h5py 
import ConfigParser
import numpy as np
import iopro
import sys

def csv2hdf5():
    # Load the configuration data. Use this for input.
    with open('./f0postProcess.cfg') as fp:
        config = ConfigParser.RawConfigParser()
        config.readfp(fp)
        files = config.get('sidechain','csvfiles')
        sc_h5file = config.get('sidechain','h5file')
        time_h5tag = config.get('sidechain','time_h5tag')
        Nsites = config.getint('sidechain','Nsites')
    
    # Load all of the CSV file paths into an array
    with open(files) as f:
        csvfiles = f.readlines()
    print("Loading data from {} csv files...".format(len(csvfiles)))
    
    # Load each CSV sequentially into an HDF file that is created on the first step
    first_file=True
    for i,file in enumerate(csvfiles):
        print "File #{}: {}".format(i, file.strip())
        x = iopro.loadtxt(file.strip(),delimiter=',')
        assert(format(x[0::Nsites].shape ==  x[(Nsites-1)::Nsites].shape))
        assert(len(x) % Nsites == 0)
        print "\tX total shape: {}".format((x.shape[0]/Nsites, Nsites, len(x[0,:])))
        Ntimes = len(x) / Nsites

        # Create the HDF file if we are looping through for the first time using the sizes from the csv files
        if first_file:
            Ncoarse = len(x[0])
            print "\tCreating new dataset, loaded file with Ncoarse = {} and Ntimes = {}".format(Ncoarse, Ntimes)
            conv.dEhdf5_init(sc_h5file, time_h5tag, 'w', Ncoarse=Ncoarse, dt_ps=.005)
            first_file=False

        # 
        h5_out =  h5py.File(sc_h5file,'a')
        ds  = h5_out[time_h5tag]
        oldlen = ds.shape[0]
        newlen = oldlen + Ntimes
        ds.resize(newlen, axis=0)
        try:
            for i in range(Nsites):
                ds[oldlen:newlen,i,:] = x[i::Nsites,:]
        except:
            print "Write failed..."
            h5_out.close()
            print("Error: {}".format(sys.exc_info()[0]))
            raise
    
def ComputeCorrelation():
    with open('./f0postProcess.cfg') as fp:
        config = ConfigParser.RawConfigParser()
        config.readfp(fp)
        sc_h5file = config.get('sidechain','h5file')
        time_h5tag = config.get('sidechain','time_h5tag')
        h5stats= config.get('sidechain','h5stats')
        h5corrtag = config.get('sidechain','h5crtag')
        h5eavtag = config.get('sidechain','h5eavtag')
    
    # Perform the computation
    sc_file = h5py.File(sc_h5file)
    try:
        sc_ds_tij = sc_file[time_h5tag]
        corr_iab, Eavg_ia = sc.AvgAndCorrelateSidechains(sc_ds_tij)
        sc_file.close()
    except:
        sc_file.close()
        print("Error: {}".format(sys.exc_info()[0]))
        raise

    # Store the computation
    corr_file = h5py.File(h5stats, 'w')
    try:
        corr_ds = corr_file.create_dataset(h5corrtag, data=corr_iab)
        Eavg_ds = corr_file.create_dataset(h5eavtag, data=Eavg_ia)
        corr_file.close()
    except:
        corr_file.close()
        print("Error: {}".format(sys.exc_info()[0]))
        raise
        

def ComputePCA():
    with open('./f0postProcess.cfg') as fp:
        config = ConfigParser.RawConfigParser()
        config.readfp(fp)
        sc_h5file = config.get('sidechain','h5file')
        time_h5tag = config.get('sidechain','time_h5tag')
        h5stats= config.get('sidechain','h5stats')
        h5eavtag = config.get('sidechain','h5eavtag')
        h5corrtag = config.get('sidechain','h5crtag')
        pca_h5file = config.get('sidechain','pcafile')

    sc_file  = h5py.File(sc_h5file)
    print sc_file.keys()
    sc_ds    = sc_file[time_h5tag]
    stat_file  = h5py.File(h5stats)

    print "Loading covariance and averages '{},{}' from hdf5 file {}...".format(h5corrtag,h5eavtag,h5stats)
    corr   = stat_file[h5corrtag]
    Eav_ij = stat_file[h5eavtag]

    print "Computing Modes..."
    eigval_in, eigvec_inj, impact_in = sc.ComputeModes(corr)
    print "Eigenvector dimension: {}".format(eigvec_inj.shape)

    conv.ApplyPCA_hdf5(sc_ds, Eav_ij, eigvec_inj, pca_h5file, time_h5tag, site=0, overwrite=True)

    sc_file.close()
    stat_file.close()

    

if __name__ == '__main__':
    csv2hdf5()
    ComputeCorrelation()
    ComputePCA()
