import depca.convert as conv
import depca.sidechain_corr as sc
import h5py 
import ConfigParser

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
    main()
