import numpy as np
import scipy.linalg as LA
import matplotlib.pylab as plt
import h5py
from copy import deepcopy
from f1SidechainCorr import h5crtag, h5eavtag
import ConfigParser
import argparse


def PlotLogSpectrum(site, vi, impact_i, floor = 1E-15, plottype = 'display', fname="spect"):
    assert(plottype != None)
    v = np.array([abs((x * imp)) + floor for x, imp in zip(vi[site], impact_i[site])])
    v = v[::-1]
    v = v[0:-3]
    #floor /= v[0]
    #v /= v[0]
    v = np.sqrt(v)
    vlog = np.log10(v)
    plt.plot(vlog, 'ro')
    #plt.plot((math.log(floor,10),) * len(vi[i]))
    DisplayPlots(plottype, fname)

def Plot2DHist(x, y, xylabels=["",""], plottype = 'display', fname = "plot2d"):
    assert len(xylabels) == 2
    xwidth = (max(x) - min(x)) / 2. * 1.10
    xmid   = (max(x) + min(x)) / 2.
    ywidth = (max(y) - min(y)) / 2. * 1.10
    ymid   = (max(y) + min(y)) / 2.
    xedges = np.linspace(xmid - xwidth, xmid + xwidth, 100)
    yedges = np.linspace(ymid - ywidth, ymid + ywidth, 100)
    
    H, xedges, yedges = np.histogram2d(x, y, bins=(xedges,yedges))
    im = plt.imshow(H, interpolation='nearest', origin='low', extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
    plt.xlabel = xylabels[0]
    plt.ylabel = xylabels[1]
    DisplayPlots(plottype, fname)

def Plot1DHist(entries, residual, plottype = 'display', fname = "plot2d", legend= [], free_energy=True, parabola=False):
    plots = []
    print entries.shape
    for i, mode in enumerate(entries):
        print "Plotting mode histogram..."
        xwidth = (max(mode) - min(mode)) / 2. * 1.10
        xmid   = (max(mode) + min(mode)) / 2.
        xedges = np.linspace(xmid - xwidth, xmid + xwidth, 100)
        H, xedges = np.histogram(mode, bins=xedges)
        xcenter = (0.5 *(xedges[1:] + xedges[:-1]))
        if free_energy:
            H = -np.log(H)
            H -= min(H) - float(i)
            p, = plt.plot(xcenter , H)
        else:
            p, = plt.plot(H, xcenter)
        plots.append(p)

    mode = residual
    xwidth = (max(mode) - min(mode)) / 2. * 1.10
    xmid   = (max(mode) + min(mode)) / 2.
    xedges = np.linspace(xmid - xwidth, xmid + xwidth, 100)
    H, xedges = np.histogram(mode, bins=xedges)
    xcenter = (0.5 *(xedges[1:] + xedges[:-1]))
    if parabola and free_energy:
        sigma_sq = np.sum( np.square(xcenter) * H / float(np.sum(H)) )
        plt.plot(xcenter, (.5 * xcenter**2 / sigma_sq) - 1)
    if free_energy:
        H = -np.log(H)
        H -= min(H) + 1
        p, = plt.plot(xcenter, H)
    else:
        p, = plt.plot(H,xcenter)
    plots.append(p)
    print plots
    print legend
    plt.legend(plots,legend)
    DisplayPlots(plottype, fname)

def PlotTimeseries(dDEmodes_nt, dDEresidual, N = None, do_sum = False, plottype='display', legend=[], fname="timeseries"):
    timeplottag = ""
    plots = []
    if N:
        spacer = dDEmodes_nt.shape[1]/N
    else:
        spacer = 1
    
    for timeseries_mode in dDEmodes_nt:
        p1, = plt.plot(timeseries_mode[::spacer])
        plots.append(p1)
    p2, = plt.plot(dDEresidual[::spacer])
    plots.append(p2)
    
    print plots
    plt.legend(plots,legend)
    DisplayPlots(plottype, fname)


def DisplayPlots(plottype, fname):
    if 'png' in plottype:
        plt.savefig("{}.png".format(fname))
    if 'pdf' in plottype:
        plt.savefig("{}.pdf".format(fname))
    if 'display' in plottype:
        plt.show()
    plt.clf()


def ComputeModes(corr, cutoffFactor = 1E-4):
	"""
	\input   corr	- a lenth Ni array of two-dimensional correlation matrix, dimension M

	\output  vi     - Ni x N eigenvalues
	         wi	- Ni x N eigenvectors of length N 
		impact_i- Ni x N values of the scale factors to conserve energy

	All outputs are sorted by vi*impact_i, with pairings preserved
	All eigenvectors are selected to have positive sum of components, and a factor of -1 is applied to those which do not.

	"""
	#print "Computing Eigenvalues..."
	vi = []
	wi = []
	Ni = len(corr)

	for i in xrange(Ni):
		v,w = LA.eigh(corr[i])
		w = w.transpose()

		valmax = max(v)
		cutpos = 0
		vi.append(np.array(v[cutpos:]))
		wi.append(np.array(w[cutpos:,:]))
	vi = np.array(vi)
	wi = np.array(wi)
	
        # NOTE: impact will take the same sign as the eigenvector, which is allowed since eigenvectors are still e-vecs under scaling.
	# impact is the mode-specific weighting factor to conserve energy under basis rotation
	impact_i = []
	for i in xrange(Ni):
		impact = []
		for v,wn in zip(vi[i], wi[i]):
			#wn /= LA.norm(wn)
			n_factor = sum(wn)
			if (v < 0 and n_factor > 0) or (v > 0 and n_factor < 0):
				wn *= -1
				n_factor *= -1
			#print n_factor
			impact.append(n_factor)
		impact_i.append(np.array(impact))
	impact_i = np.array(impact_i)
	# --------------------------SORT-------------------------------

	for i in xrange(Ni):
		wi_copy = deepcopy(wi[i])
		tup = zip(vi[i], impact_i[i], wi_copy)
		tup.sort(key=lambda x: x[0] * x[1])
		for j in xrange(len(tup)):
			vi[i,j] = tup[j][0]
			impact_i[i,j] = tup[j][1]
			wi[i,j] = tup[j][2]

	return vi, wi, impact_i

def DeltaModeTracker(E_t_ij, E_avg_ij, modes_inj, site, modes_requested=[], Nframes=None, offset=0):
	"""
	DeltaModeTracker: Takes E(t) np.ndarray and an iterable of modes, and plots the deviation of dE(t) [after taking the average for each site]
	                  as contributed from the modes, as compared to the actual dE(t)

	\input 	E_t_ij 		- a collection of Eij as a function of time
		modes_inj 	- a collection of n normalized modes (in j coordinates) for each site i -- rotation matrix between bases

	\output dDE_t	        - a list, for (site,mode) pairs (i,n) requested, of 1d arrays tracking those modes in time
	"""

	print "Running DeltaModeTracker..."
	Ntimes  = E_t_ij.shape[0]
	Nsites  = E_t_ij.shape[1]
	Ncoarse = E_t_ij.shape[2]

        # Flip the modes so that they are in order from largest to smallest and zero-based 

	E_avg_ij = np.zeros( (Nsites, Ncoarse) )
	print "Ntimes: ", Ntimes
	print "Nsites: ", Nsites
	print "Ncoarse: ", Ncoarse
        print "Requested modes (greatest-to-largest, one-based-indexing):", modes_requested
        modes_requested = [Ncoarse - mode for mode in modes_requested]
	
	# First, check that the modes are normalized, then weight the modes
	print "Checking normalization and applying weights..."
        modeweight_inj = np.zeros(modes_inj.shape)
        modeweight_inj[:] = modes_inj[:]
	for m_nj in modeweight_inj:
            for m_j in m_nj:
                if not abs(LA.norm(m_j) - 1) < 1E-6:
                    raise RuntimeError("Bad normalization encountered, {}".format(LA.norm(m_j)-1))
                # Why are my modes squaring??? In the mode matrix???
                m_j *= sum(m_j)
	
        if Nframes:
            outLen = min(Ntimes-offset,Nframes)
        else:
            outLen = Ntimes
        print outLen
    
        return_modes_nt = np.zeros((len(modes_requested), outLen))
	#Then, compute mode timeseries:
	print "Running vectorized mode tracker computation (O(T*num(modes_requested))) using the total mean...,"
        #dE_t_j = np.zeros((outLen,E_t_ij.shape[2]))
        dE_t_j = E_t_ij[offset:offset+outLen,site,:] - np.mean(E_t_ij[:,site,:], axis=0)
        print "dE matrix completed...",
	for i,n in enumerate(modes_requested):
            print "mode {} completed...".format(n),
            return_modes_nt[i,:] = np.inner(dE_t_j[0:outLen,:], modeweight_inj[site,n,:])
        print 'done'


	print "Running residual dDE(t) computation (O(T))..."
        dDEresidual_t = (np.sum(dE_t_j[:outLen,:], axis=1) - np.sum(return_modes_nt, axis=0))

	print "Done."

	return return_modes_nt, dDEresidual_t



#--------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Compute the eigenvalues of the correlation matrix produced by SidechainCorr, then plot the timeseries for the modes selected. Leaves the database unmodified.")
    parser.add_argument("site", type=int, help="Site that the data is requested for, use 1-based indexing. No error checking.")
    parser.add_argument("-dEtmodes", type=int, nargs='+', action='append',     help="(requires plot dEt) A collection of all modes to include in the timeseries, using zero-based indexing.")
    parser.add_argument("-modes2dhist", type=int, nargs=2, action='append',        help="(requires plot2dhist) Selecting the mode-pair for 2d histogram. Repeat this option for multiple plots.")
    parser.add_argument("-modes1dhist", type=int, nargs='*', action='append',       help="(requires plot1dhist) Select modes to histogram together. Modes in the same option will be plotted together. Repeat this option for multiple plots.")
    parser.add_argument("-outfnamebase", default="plot", help="Filename base for output from plotspectrum and other outputs")
    parser.add_argument("-plotspectrum", action='store_true', help="Set to plot PCA spectrum")
    parser.add_argument("-parabola", action='store_true', help="Compare 1-D histograms to their mean/variance parabola.")
    parser.add_argument("-savemode", default=['display'], nargs='+', choices=['pdf','png','display'], help="Set format for file output")
    parser.add_argument("-Nframes", type=int, help="Number of frames to include in calculations")
    parser.add_argument("-offset", type=int, default=0, help="Number of frames to skip before beginning computation")


    args = parser.parse_args()
    print args
    args.site -= 1
    print args.savemode


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
        vi, wi, impact_i = ComputeModes(corr)
        print "\tModes Computed!"
    
        if args.plotspectrum:
            print "Plotting spectrum..."
            fname = "{}_spect{}".format(args.outfnamebase,args.site+1)
            PlotLogSpectrum(args.site, vi, impact_i, plottype=args.savemode, fname=fname)
        
        if args.modes2dhist:
            for modepair in args.modes2dhist:
                print "Plotting 2D histogram for site {}...".format(args.site+1)
                xylabels = ["Mode {}".format(modepair[0]), "Mode {}".format(modepair[1])]
                fname = "{}2D_s{}_{}v{}".format(args.outfnamebase, args.site+1, modepair[0], modepair[1])
                dDE_nt,_ = DeltaModeTracker(E_t_ij, Eav_ij, wi, args.site, modepair, Nframes=args.Nframes, offset=args.offset)
                Plot2DHist(dDE_nt[0], dDE_nt[1], xylabels = xylabels, plottype = args.savemode, fname=fname)
        
        if args.dEtmodes:
            for modeset in args.dEtmodes:
                legend = ["Mode {}".format(mode) for mode in modeset]
                legend.append("Residual")
                print "LEGEND: ", legend
                fname = "{}dEt_s{}".format(args.outfnamebase, args.site+1)
                dDE_nt, residual_t = DeltaModeTracker(E_t_ij, Eav_ij, wi, args.site, modeset, Nframes=args.Nframes, offset=args.offset)
                PlotTimeseries(dDE_nt, residual_t, plottype = args.savemode, legend=legend, fname=fname)

        if args.modes1dhist:
            for modeset in args.modes1dhist:
                legend = ["Mode {}".format(mode) for mode in modeset]
                legend.append("Residual")
                print "LEGEND: ", legend
                fname = "{}dEt_s{}".format(args.outfnamebase, args.site+1)
                dDE_nt, residual_t = DeltaModeTracker(E_t_ij, Eav_ij, wi, args.site, modeset, Nframes=args.Nframes, offset=args.offset)
                Plot1DHist(dDE_nt, residual_t, legend=legend, plottype=args.savemode, fname=fname, parabola=args.parabola)



if __name__ == "__main__":
	main()
