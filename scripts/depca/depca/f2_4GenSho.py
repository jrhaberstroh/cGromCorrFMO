import numpy as np
import scipy.linalg as LA
import matplotlib.pylab as plt
import h5py
from copy import deepcopy
from f1SidechainCorr import h5crtag, h5eavtag
import ConfigParser
import argparse
from f2_1TrackModes import ComputeModes, DeltaModeTracker, ComputeCt_FFT, DisplayPlots


def SHOdata(dDE_t, dt, fname=None):
    Xsq = np.mean(np.square(dDE_t))
    Ct = ComputeCt_FFT(dDE_t, dDE_t, len(dDE_t)/2 )
    Ct /= Ct[0]
    # Compute the first crossing, or use tf = inf if no crossings
    crossings = np.nonzero(Ct[1:]*Ct[:-1] <= 0.0)
    if len(crossings) > 0:
        tf = crossings[0][0] + 1
    else:
        tf = len(dDE_t)
    # Compute the time constant as the integral under the curve
    tau = np.sum(Ct[0:tf])*dt
    print "{},{} ".format(tau,Xsq)
    plt.plot(np.arange(tf)*dt, Ct[0:tf])
    if fname:
        DisplayPlots(['png'],fname)
    return tau, Xsq

#--------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute the eigenvalues of the correlation matrix produced by SidechainCorr, then plot the timeseries for the modes selected. Leaves the database unmodified.")
    parser.add_argument("modes", type=int, nargs='*', help="modes to generate SHOs for, from [1, nModes] in [max,min]")
    parser.add_argument("-site", type=int, default=1, help="Site that the data is requested for, use 1-based indexing. No error checking.")
    parser.add_argument("-Nframes", type=int, help="Number of frames to include in calculations")
    parser.add_argument("-offset", type=int, default=0, help="Number of frames to skip before beginning computation")

    args = parser.parse_args()
    print args
    args.site -= 1

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
        dt = E_t_ij.attrs["dt"]
    
        print "Computing Modes..."
        vi, wi, impact_i = ComputeModes(corr)
        print "\tModes Computed!"

        if len(args.modes) == 0:
            args.modes = np.arange(wi.shape[-1]) + 1
            print "No modes provided, computing all modes"
            print args.modes
    
        dDE_nt,_ = DeltaModeTracker(E_t_ij, Eav_ij, wi, args.site, args.modes, Nframes=args.Nframes, offset=args.offset)
        for n in xrange(dDE_nt.shape[0]):
            fname = "../../data/images/Ct_site{}_{}.png".format(args.site+1, n+1)
            SHOdata(dDE_nt[n,:], dt, fname)


if __name__=="__main__":
    main()
