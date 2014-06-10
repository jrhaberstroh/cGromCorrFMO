import numpy as np
import numpy.linalg as LA
import h5py
import matplotlib.pyplot as plt
import argparse
import ConfigParser
from f2_1TrackModes import ComputeModes, DeltaModeTracker
from f1SidechainCorr import h5crtag, h5eavtag


def main():
    parser = argparse.ArgumentParser(description="Analyze anharmonicity of dynamics through various lenses")
    parser.add_argument("-site", type=int, default=1, help="Site that the data is requested for, use 1-based indexing. No error checking.")
    parser.add_argument('-computesimga', action='store_true', help="Compute and plot the percentage of the residual sigma that comes from the modes")
    args = parser.parse_args()


    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')
    h5file = config.get('sidechain','h5file')
    timetag = config.get('sidechain','time_h5tag')
    corrtag = config.get('sidechain','corr_h5tag')


    with h5py.File(h5file,'r') as f:
        print "Loading dset '{}' from hdf5 file {}...".format(timetag,h5file)
        E_t_ij = f[timetag]
        corr   = f[corrtag+h5crtag]
        Eav_ij = f[corrtag+h5eavtag]
    
        print "Computing Modes..."
        wi, vi, impact_i = ComputeModes(corr, sort=True)
        print "\tModes Computed!"
    
        evs = vi[args.site-1,:,:] 
        strs = wi[args.site-1,:] * np.square(impact_i[args.site-1,:])
        plt.plot(np.sort(strs)[::-1], linewidth=2)
        plt.plot(wi[args.site-1,:][::-1], linewidth=2)
        #plt.semilogy(np.sort(strs), linewidth=2)
        #plt.semilogy(wi[args.site-1,:], linewidth=2)
        plt.show()
        plt.imshow(np.inner(evs.T, evs.T), interpolation='none', cmap='gray')
        plt.colorbar()
        plt.show()

        s = args.site-1
        
        sigma    = np.zeros(wi[s,:].shape)
        sigmares = np.zeros(wi[s,:].shape)
        
        #print "Computing total variance...",
        #sigmatot = np.var(np.sum(E_t_ij[:,s,:], axis=1))
        #print "done."

        E_t_mode = np.zeros(E_t_ij.shape[0])
        E_t_res  = np.zeros(E_t_ij.shape[0])
        for mu in xrange(sigma.shape[0]):
            print mu
            E_t_mode[:] = np.inner(E_t_ij[:,s,:], vi[s,:,mu]) * impact_i[s,mu]
            E_t_res += E_t_mode
            sigma[mu]    = np.var(E_t_mode)
            sigmares[mu] = np.var(E_t_res)


            
        




if __name__=="__main__":
    main()
