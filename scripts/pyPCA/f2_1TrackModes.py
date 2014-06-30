import numpy as np
import scipy.linalg as LA
import matplotlib.pylab as plt
import h5py
from copy import deepcopy
from f1SidechainCorr import h5crtag, h5eavtag
import ConfigParser
import argparse


def DisplayPlots(plottype, fname):
    if 'png' in plottype:
        plt.savefig("{}.png".format(fname))
    if 'pdf' in plottype:
        plt.savefig("{}.pdf".format(fname))
    if 'display' in plottype:
        plt.show()
    plt.clf()


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

    def Plot1DHistMode(mode, offset, plots, free_energy=True, parabola=False):
        print "Plotting mode histogram..."
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
            H -= min(H) - float(offset)
            p, = plt.plot(xcenter , H)
        else:
            p, = plt.plot(H, xcenter)
        plots.append(p)


    for i, mode in enumerate(entries):
        Plot1DHistMode(mode, offset=i, plots=plots, free_energy=free_energy, parabola=False)

    Plot1DHistMode(residual, offset=-1, plots=plots, free_energy=free_energy, parabola=False)

    print plots
    print legend
    plt.legend(plots,legend)
    plt.xlabel('dE / cm^1')
    plt.ylabel('-ln(p)')
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


def PlotCt(dDEmodes_nt, dDEresidual, N = None, totalCt=False, plottype='display', legend=[], fname="timeseries", dt = 1.):
    timeplottag = ""
    plots = []
    ctlen = dDEmodes_nt.shape[1] / 2
    t = np.linspace(0, (ctlen-1)*dt, ctlen)
    
    def PlotCtMode(mode_t, ctlen):
        Ct = ComputeCt_FFT(mode_t,mode_t, ctlen)
        Ct/= Ct[0]
        p, = plt.plot(t, Ct)
        plots.append(p)

    for mode_t in dDEmodes_nt:
        PlotCtMode(mode_t, ctlen)
    PlotCtMode(dDEresidual, ctlen)
    if totalCt:
        PlotCtMode(dDEresidual + np.sum(dDEmodes_nt,axis=0), ctlen)
    
    ax = plt.gca()
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    plt.legend(plots,legend)
    plt.xlabel = "Time, ps"
    plt.ylabel = "Correlation, scaled"
    DisplayPlots(plottype, fname)


def ComputeModes(corr, sort=True):
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
        if sort:
	    for i in xrange(Ni):
		    wi_copy = deepcopy(wi[i])
		    tup = zip(vi[i], impact_i[i], wi_copy)
		    tup.sort(key=lambda x: x[0] * x[1] * x[1])
		    for j in xrange(len(tup)):
			    vi[i,j] = tup[j][0]
			    impact_i[i,j] = tup[j][1]
			    wi[i,j] = tup[j][2]

	return vi, wi, impact_i

def ComputeCt_FFT(f1_t, f2_t, corr_len=None):
    if len(f1_t) != len(f2_t):
        raise TypeError("ComputeCt trying to correlate functions of different sizes: {} and {}".format(len(f1_t), len(f2_t)))
    if not corr_len:
        corr_len = len(f1_t)
    numavg  = len(f1_t) - np.arange(corr_len)
    mycorr = np.real(np.fft.ifft( np.fft.fft(f2_t) * np.conj(np.fft.fft(f1_t)) )[0:corr_len])/numavg
    return mycorr

def DeltaModeTracker(E_t_ij, E_avg_ij, modes_inj, site, modes_requested=[], Nframes=None, offset=0):
	"""
	DeltaModeTracker: Takes E(t) np.ndarray and an iterable of modes, and plots the deviation of dE(t) [after taking the average for each site]
	                  as contributed from the modes, as compared to the actual dE(t)

	\input 	E_t_ij 		- a collection of Eij as a function of time
                E_avg_ij        - open HDF5 dataset of average energy on mode j of site i
		modes_inj 	- a collection of n normalized modes (in j coordinates) for each site i -- rotation matrix between bases
                modes_requested - (opt) a 1-based list of modes to compute timeseries for (e.g. modes_requested=[1,2,3])

	\output dDE_t	        - a list, for (site,mode) pairs (i,n) requested, of 1d arrays tracking those modes in time
	"""

	print "Running DeltaModeTracker..."
	Ntimes  = E_t_ij.shape[0]
	Nsites  = E_t_ij.shape[1]
	Ncoarse = E_t_ij.shape[2]

        # Flip the modes so that they are in order from largest to smallest and zero-based 
	print "Ntimes: ", Ntimes
	print "Nsites: ", Nsites
	print "Ncoarse: ", Ncoarse
        print "Requested modes (greatest-to-largest, one-based-indexing):", modes_requested
        modes_requested = [Ncoarse - mode for mode in modes_requested]

        # Compute t0 and tf from Ntimes and offset
        if Nframes:
            outLen = min(Ntimes-offset,Nframes)
        else:
            outLen = Ntimes
        print outLen
        t0 = offset
        tf = offset+outLen

	# First, check that the modes are normalized, then weight the modes
	print "Checking normalization and applying weights..."
        modeweight_inj = np.zeros(modes_inj.shape)
        modeweight_inj[:] = modes_inj[:]
	for m_nj in modeweight_inj:
            for m_j in m_nj:
                if not abs(LA.norm(m_j) - 1) < 1E-5:
                    raise RuntimeError("Bad normalization encountered, {}".format(LA.norm(m_j)-1))
                # Why are my modes squaring??? In the mode matrix???
                m_j *= sum(m_j)
	
    
	#Then, compute mode timeseries:
	print "Computing dE for site {} (O(T*num(modes_requested))) using the total mean...".format(site+1)

        GB = float(1E9)
        RAM_GB = 2
        RAM_nfloat = RAM_GB * 1E9 / 8
        RAM_ntimes = RAM_nfloat / Ncoarse
        RAM_nchunk = int( np.ceil((tf - t0) / float(RAM_ntimes)) )
        RAM_time_per_chunk = (tf - t0) / RAM_nchunk
        print "Number of chunks needed: {} of {}GB each".format(RAM_nchunk, RAM_time_per_chunk * Ncoarse * 8 / GB)
        RAM_return_times = (RAM_nchunk*RAM_time_per_chunk)
        RAM_return_tot = (RAM_nchunk*RAM_time_per_chunk) * 8 * (len(modes_requested) + 1)
        print "RAM needed for return data: {}GB".format(RAM_return_tot / GB)


        dE_t_j = np.zeros((RAM_time_per_chunk,Ncoarse))
        return_modes_nt = np.zeros((len(modes_requested), RAM_return_times))
        dDEresidual_t = np.zeros((RAM_return_times))
        #TODO: Implement chunking to compute large datasets
        for chunk_num in xrange(RAM_nchunk):
            print "Chunk {}:".format(chunk_num+1)
            # Each chunk reads [t0_chunk, tf_chunk)
            t0_chunk = (  chunk_num   * RAM_time_per_chunk) + t0
            tf_chunk = ((chunk_num+1) * RAM_time_per_chunk) + t0
            t0_return = t0_chunk - t0
            tf_return = tf_chunk - t0

            # Build dE for chunk
            print "Computing chunk dE complete...",
            dE_t_j[:,:] = E_t_ij[ t0_chunk:tf_chunk, site,:] - E_avg_ij[site,:]
            print "chunk dE complete..."

	    for i,n in enumerate(modes_requested):
                print "Computing mode {}...".format(n+1),
                
                rotated_modes = np.inner(dE_t_j[:,:], modeweight_inj[site,n,:])
                return_modes_nt[i,t0_return:tf_return] = rotated_modes

                print "mode {} computed...".format(n+1)

	    print "Running residual dDE(t) computation...",
            dDEresidual_t = (np.sum(dE_t_j[:,:], axis=1) - np.sum(return_modes_nt[:,t0_return:tf_return], axis=0))
            print 'Chunk {}  done'.format(chunk_num+1)

	print "Done."
	return return_modes_nt, dDEresidual_t



#--------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Compute the eigenvalues of the correlation matrix produced by SidechainCorr, then plot the timeseries for the modes selected. Leaves the database unmodified.")
    parser.add_argument("-site", type=int, default=1,                                                   help="Site that the data is requested for, use 1-based indexing. No error checking.")
    parser.add_argument("-dEtmodes", type=int, nargs='+', action='append',                              help="(requires plot dEt) A collection of all modes to include in the timeseries, using zero-based indexing.")
    parser.add_argument("-Nframes", type=int,                                                           help="Number of frames to include in calculations")
    parser.add_argument("-offset", type=int, default=0,                                                 help="Number of frames to skip before beginning computation")
    parser.add_argument("-plotspectrum", action='store_true',                                           help="Set to plot PCA spectrum")
    parser.add_argument("-hist2d", type=int, nargs=2, action='append',                                  help="A list of mode-pairs for 2d histogram. Repeat option for multiple plots.")
    parser.add_argument("-hist2d_opt", type=int, nargs=2, action='append',                              help="Options for 2D hist. None available.")
    parser.add_argument("-hist1d", type=int, nargs='*', action='append',                                help="A list of modes to histogram together. Modes will be plotted with residual energy distribution. Repeat option for multiple plots.")
    parser.add_argument("-hist1d_opt", type=str, nargs='*', default=[], choices=['parabola'],           help="Options for 1D hist. Parabola plots a parabola with the residual dE.")
    parser.add_argument("-modesCt", type=int, nargs='*', action='append',                               help="A list of modes to compute time-correlations for. Modes will be plotted with residual correlation. Repeat option for multiple plots.")
    parser.add_argument("-outfnamebase", default="plot",                                                help="Filename base for output from plotspectrum and other outputs")
    parser.add_argument("-savemode", default=['display'], nargs='+', choices=['pdf','png','display'],   help="Set format for file output")
    parser.add_argument("-evmtx", action='store_true',                                                  help="Analyze the eigenvalue matrix")


    args = parser.parse_args()
    print args
    args.site -= 1
    print args.savemode


    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')
    h5file = config.get('sidechain','h5file')
    h5stats= config.get('sidechain','h5stats')
    timetag = config.get('sidechain','time_h5tag')
    corrtag = config.get('sidechain','corr_h5tag')
    
   
    with h5py.File(h5file,'r') as f_time:
        if (h5file != h5stats):
            f_corr =  h5py.File(h5stats,'r')
        else:
            f_corr = f_time
        print "Loading timeseries '{}' from hdf5 file {}...".format(timetag,h5file)
        E_t_ij = f_time[timetag]
        print "Loading covariance and averages '{},{}' from hdf5 file {}...".format(corrtag+h5crtag,corrtag+h5eavtag,h5stats)
        corr   = f_corr[corrtag+h5crtag]
        Eav_ij = f_corr[corrtag+h5eavtag]
    
        print "Computing Modes..."
        vi, wi, impact_i = ComputeModes(corr)
        print "\tModes Computed!"
    
        if args.plotspectrum:
            print "Plotting spectrum..."
            fname = "{}_spect{}".format(args.outfnamebase,args.site+1)
            PlotLogSpectrum(args.site, vi, impact_i, plottype=args.savemode, fname=fname)

        if args.evmtx:
            raise NotImplementedError("No analysis tools for the eigenvector matrix is available yet")
        
        print "1D HIST MODES: ", args.hist1d
        if args.hist1d:
            for modeset in args.hist1d:
                legend = ["Mode {}".format(mode) for mode in modeset]
                legend.append("Residual")
                print "LEGEND: ", legend
                fname = "{}_s{}_1d{}".format(args.outfnamebase, args.site+1, '-'.join([str(j) for j in modeset]))
                dDE_nt, residual_t = DeltaModeTracker(E_t_ij, Eav_ij, wi, args.site, modeset, Nframes=args.Nframes, offset=args.offset)
                Plot1DHist(dDE_nt, residual_t, legend=legend, plottype=args.savemode, fname=fname, parabola='parabola' in args.hist1d_opt)

        if args.hist2d:
            for modepair in args.hist2d:
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


        if args.modesCt:
            for modeset in args.modesCt:
                legend = ["Mode {}".format(mode) for mode in modeset]
                legend.append("Residual")
                legend.append("Total")
                print "LEGEND: ", legend
                fname = "{}_s{}_Ct{}".format(args.outfnamebase, args.site+1, '-'.join([str(j) for j in modeset]))
                dDE_nt, residual_t = DeltaModeTracker(E_t_ij, Eav_ij, wi, args.site, modeset, Nframes=args.Nframes, offset=args.offset)
                print "Modes retrieved, computing Ct..."
                PlotCt(dDE_nt, residual_t, legend=legend, totalCt = True, plottype=args.savemode, fname=fname, dt = E_t_ij.attrs["dt"])

        f_corr.close()


if __name__ == "__main__":
	main()
