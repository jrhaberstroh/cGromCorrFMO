import depca.convert as conv
import depca.sidechain_corr as sc
import h5py 
import ConfigParser
import numpy as np
import iopro

def new_sidechain():
    with open('./f0postProcess.cfg') as fp:
        config = ConfigParser.RawConfigParser()
        config.readfp(fp)
        files = config.get('sidechain','csvfiles')
        sc_h5file = config.get('sidechain','h5file')
        time_h5tag = config.get('sidechain','time_h5tag')
    
    with open(files) as f:
        csvfiles = f.readlines()
    print("Loading data from {} csv files...".format(len(csvfiles)))

    first_file=True
    Nsites = 7
    for i,file in enumerate(csvfiles):
        print "File #{}: {}".format(i, file.strip())
        adapter = iopro.text_adapter(file.strip(),parser='csv', field_names=False)
        x = adapter[:][:]
        print "SHAPE OF LOADED DATA: {}".format(x.shape)
        print "SHAPE OF FIRST INDEX: {}".format(len(x[0]))
        Ntimes = len(x) / Nsites
        if first_file:
            Ncoarse = len(x[0])
            print "Creating new dataset, loaded file with Ncoarse = {} and Ntimes = {}".format(Ncoarse, Ntimes)
            conv.dEhdf5_init(sc_h5file, time_h5tag, 'w', Ncoarse=Ncoarse, dt_ps=.005)
            first_file=False

        h5_out =  h5py.File(sc_h5file,'a')
        ds  = h5_out[time_h5tag]
        oldlen = ds.shape[0]
        newlen = oldlen + Ntimes
        ds.resize(newlen, axis=0)
        for i in range(Nsites):
            ds[oldlen:newlen,i,:] = x[i::7,:]
        h5_out.close()




    


def main():
    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')

    sc_h5file = config.get('sidechain','h5file')
    pca_h5file = config.get('sidechain','pcafile')
    h5stats= config.get('sidechain','h5stats')
    corrtag = config.get('sidechain','corr_h5tag')
    h5crtag = config.get('sidechain','h5crtag')
    time_h5tag = config.get('sidechain','time_h5tag')
    h5eavtag = config.get('sidechain','h5eavtag')
    config.close()

    sc_file  = h5py.File(sc_h5file)
    sc_ds    = sc_file[time_h5tag]
    stat_file  = h5py.File(h5stats)

    print "Loading covariance and averages '{},{}' from hdf5 file {}...".format(corrtag+h5crtag,corrtag+h5eavtag,h5stats)
    corr   = stat_file[corrtag+h5crtag]
    Eav_ij = stat_file[corrtag+h5eavtag]

    print "Computing Modes..."
    eigval_in, eigvec_inj, impact_in = sc.ComputeModes(corr)
    print "Eigenvector dimension: {}".format(eigvec_inj.shape)

    conv.ApplyPCA_hdf5(sc_ds, Eav_ij, eigvec_inj, pca_h5file, time_h5tag, site=0, overwrite=True)

    sc_file.close()
    stat_file.close()

    

if __name__ == '__main__':
    new_sidechain()
    #main()
