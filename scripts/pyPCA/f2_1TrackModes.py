import numpy as np
import cPickle
import scipy.linalg as LA
import matplotlib
import matplotlib.pylab as plt
import math
import sys
import h5py
from copy import deepcopy
from f1SidechainCorr import h5tstag, h5crtag
import ConfigParser
import argparse

# -----------------------------OLD CODE :(---------------------------------------

def dEModesSeparate(dEmodes_i_tn, dEtot_t_i, n_modes, i = 0):
	assert(type(dEmodes_i_tn) == np.ndarray)
	assert(n_modes <= dEmodes_i_tn.shape[2])
	Ntimes = dEmodes_i_tn.shape[1]

	timeseries_sum = np.zeros(Ntimes)

	plots = []
	legend = []
	for mode_counter in xrange(n_modes):
		mode = -mode_counter - 1
		timeseries_mode = dEmodes_i_tn[i,:,mode]

		p1, = plt.plot(timeseries_mode)
		plots.append(p1)
		legend.append("Mode "+str(mode_counter+1))


	plt.legend(plots,legend)
	plt.savefig("plot_timesep_m"+str(n_modes)+"s"+str(i+1)+".png")
	plt.savefig("plot_timesep_m"+str(n_modes)+"s"+str(i+1)+".pdf")
	#plt.show()




def ModeTracker(E_t_ij, modes_inj):
	"""
	\input 	E_t_ij 		- a collection of Eij as a function of time
		modes_inj 	- a collection of n modes (in j coordinates) for each site i

	\output Emodes_t_in	- a collection of the energetics of the n modes for site i as a function of time
	"""

	Emodes_t_in = []

	for t, Eij in enumerate(E_t_ij):
		if t%100 == 0:
			print "t =",t
		Emodes_t_in.append([])
		for i, Ej in enumerate(Eij):
			Emodes_t_in[t].append([])
			for n, Mj in enumerate(modes_inj[i]):
				Emodes_t_in[t][i].append(np.dot(Ej, Mj))

	return Emodes_t_in

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
	#print "\tComputed"
	
	
	#print "Checking normalization, computing normalized impact, and converting all impacts to positive numbers..."
	# impact will take the same sign as the eigenvector, which is allowed since eigenvectors are still e-vecs under scaling.
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

		

#		vi[i], impact_i[i], wi[i] = zip(*sorted(zip(vi[i],impact_i[i], wi[i]), key=lambda val: val[0] * val[1]))
#		#print "EIGENVALUES: ", srtd[:][0], "length = ", len(srtd[:][0])
#		#impact_i = x[1]
	#print "\tNormalization imposed"
	return vi, wi, impact_i


# -----------------------------NEW CODE!---------------------------------------

def PlotLogSpectrum(site, vi, impact_i, floor = 1E-15, plottype = None, fnamebase="plot"):
    assert(plottype != None)
    v = np.array([abs((x * imp)) + floor for x, imp in zip(vi[i], impact_i[i])])
    v = v[::-1]
    v = v[0:-3]
    #floor /= v[0]
    #v /= v[0]
    v = np.sqrt(v)
    vlog = np.log10(v)
    plt.plot(vlog, 'ro')
    #plt.plot((math.log(floor,10),) * len(vi[i]))

    if plottype == 'pdf':
        plt.savefig("{}_spect{}.pdf".format(fnamebase,site))
        plt.clf()
    if plottype == 'png':
        plt.savefig("{}_spect{}.png".format(fnamebase,site))
        plt.clf()
    if plottype == 'display':
        plt.show()

def Plot2DHist(x, y, site=None, mode1=None, mode2=None, plottype = None, fnamebase = "plot"):
    if fnamebase == "":
        fnamebase = "plot"
    xwidth = (max(x) - min(x)) / 2. * 1.10
    xmid   = (max(x) + min(x)) / 2.
    ywidth = (max(y) - min(y)) / 2. * 1.10
    ymid   = (max(y) + min(y)) / 2.
    xedges = np.linspace(xmid - xwidth, xmid + xwidth, 50)
    yedges = np.linspace(ymid - ywidth, ymid + ywidth, 50)
    
    print "Saving histogram figure..."
    H, xedges, yedges = np.histogram2d(x, y, bins=(xedges,yedges))
    im = plt.imshow(H, interpolation='nearest', origin='low', extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
    if plottype=='png':
        plt.savefig("{}_2d_s{}_{}v{}.png".format(fnamebase,site,mode1,mode2))
        plt.clf()
    if plottype=='pdf':
        plt.savefig("{}_2d_s{}_{}v{}.pdf".format(fnamebase,site,mode1,mode2))
        plt.clf()
    if plottype=='display':
        plt.show()

def PlotDeltaTimeseries(dEmodes_i_tn, dEtot_t_i, n_modes, i = 0):
	assert(type(dEmodes_i_tn) == np.ndarray)
	assert(n_modes <= dEmodes_i_tn.shape[2])
	Ntimes = dEmodes_i_tn.shape[1]

	timeseries_sum = np.zeros(Ntimes)
	for mode_counter in xrange(n_modes):
		mode = -mode_counter - 1
		timeseries_mode = dEmodes_i_tn[i,:,mode]

		timeseries_sum += timeseries_mode

	timeseries_sum -= dEtot_t_i[:,i]

	p1, = plt.plot(timeseries_sum)
	p2, = plt.plot(dEtot_t_i[:,i])
	plt.legend([p1,p2],["Difference from dEtot", "dEtot"])
	#plt.plot(dEtot_t_i[:,i])
	plt.savefig("plot_timetot_m"+str(n_modes)+"s"+str(i+1)+".png")
	plt.savefig("plot_timetot_m"+str(n_modes)+"s"+str(i+1)+".pdf")
	#plt.show()

def DeltaModeTracker(E_t_ij, modes_inj, site, modes_requested=[]):
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
	
	print "Computing time-average for each site O(T)..."
        E_avg_ij = np.mean(E_t_ij, axis=0)

    
	# First, check that the modes are normalized, then weight the modes
	print "Checking normalization and applying weights..."
	for m_nj in modes_inj:
		for m_j in m_nj:
			assert abs(LA.norm(m_j) - 1) < 1E-6
			m_j *= sum(m_j)
	
        return_modes_nt = np.zeros((len(modes_requested), Ntimes))
	#Then, compute mode timeseries:
	print "Running vectorized mode tracker computation (O(T*num(modes_requested)))..."
	for i,n in enumerate(modes_requested):
            return_modes_nt[i,:] = np.inner(E_t_ij[:,site,:], modes_inj[site,n,:]) -\
                                   np.inner(E_avg_ij[site,:], modes_inj[site,n,:])


	print "Running residual dDE(t) computation (O(T))..."
        dDEresidual_t = np.sum(E_t_ij[:,site,:], axis=1) - np.sum(E_avg_ij[site,:]) - np.sum(return_modes_nt, axis=0)

	print "Done."

	return return_modes_nt, dDEresidual_t








#--------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Compute the eigenvalues of the correlation matrix produced by SidechainCorr, then plot the timeseries for the modes selected. Leaves the database unmodified.")
    parser.add_argument("site", type=int, help="Site that the data is requested for, use 1-based indexing. No error checking.")
    parser.add_argument("--outfnamebase", help="Filename base for output from plotspectrum and other outputs")
    parser.add_argument("--plotspectrum", choices=[None,'pdf','png','display'], help="Select destination for plotted spectrum")
    parser.add_argument("--plotdEt",      choices=[None,'pdf','png','display'], help="Select destination for plotted spectrum")
    parser.add_argument("--dEtmodes", type=int, nargs='+', action='append',     help="(requires plot dEt) A collection of all modes to include in the timeseries, using zero-based indexing.")
    parser.add_argument("--plot2dhist",   choices=[None,'pdf','png','display'], help="Select destination for anharmonic scatterplot")
    parser.add_argument("--modes2d", type=int, nargs=2, action='append',        help="(requires plot2dhist) Selecting the mode-pair for 2d histogram. Repeat this option for multiple plots.")
    parser.add_argument("--plot1dhist",   choices=[None,'pdf','png','display'], help="Select destination for anharmonic scatterplot")
    parser.add_argument("--modes1d", type=int, nargs='+', action='append',       help="(requires plot1dhist) Select modes to histogram together. Modes in the same option will be plotted together. Repeat this option for multiple plots.")

    args = parser.parse_args()
    print args
    args.site -= 1

    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')
    h5file = config.get('sidechain','h5file')
    timetag = config.get('sidechain','time_h5tag')
    corrtag = config.get('sidechain','corr_h5tag')
    
    # NOTE: From old version
    #save_mode_ts = config.getboolean('modes','save_mode_ts')
    #plot_spectrum = config.getboolean('modes','plot_spectrum')
    #modes_pkl = config.get('modes','modes_pkl')
    #time_max = config.getint('modes','time_max')
    
    
    with h5py.File(h5file,'r') as f:
        print "Loading dset '{}' from hdf5 file {}...".format(timetag,h5file)
        E_t_ij = f[timetag]
        corr   = f[corrtag+h5crtag]
    
        print "Computing Modes..."
        vi, wi, impact_i = ComputeModes(corr)
        print "\tModes Computed!"
    
        if args.plotspectrum:
            print "Plotting spectrum..."
            if args.outfnamebase:
       	        PlotLogSpectrum(site, vi, impact_i, plottype=args.plotspectrum, fnamebase=args.outfnamebase)
            else:
       	        PlotLogSpectrum(site, vi, impact_i, plottype=args.plotspectrum)
        
        if args.plot2dhist:
            for modepair in args.modes2d:
                print "Plotting 2D histogram for site {}...".format(site)
                dDE_nt,_ = DeltaModeTracker(E_t_ij, wi, args.site, modepair)
                Plot2DHist(dDE_nt[0], dDE_nt[1], site=site, mode1=modepair[0], mode2=modepair[2], plottype = args.plot2dhist, fnamebase=args.outfnamebase)
               
        if args.dEtmodes:
            dDE_nt, residual_t = DeltaModeTracker(E_t_ij, wi, args.site, args.dEtmodes)

        if args.plot1dhist:
            raise NotImplementedError("plot1dhist not yet available")


    
    # NOTE: From old version
    #if not save_mode_ts:
    #    print "save_mode_ts == no, exiting TrackModes.main()"
    #    exit(101)
    #else:
    #    print "save_mode_ts = yes, saving the computed trajectories for PlotModeTimeseries..."
    #N = time_max
    #E_t_i = []
    #nread = 0
    #
    #print "Reading",N,"times of sidechain data from "+csv_filename+" to plot timeseries..."
    #E_t_i, nread = SidechainRead(csv_filename, N)
    #print "\tActually read", nread, "times from the file."
    #
    #DeltaEmodes_i_tn, DeltaEtot_t_i = np.array(DeltaModeTracker(E_t_i, wi))
    #DeltaEmodes_i_tn = np.array(DeltaEmodes_i_tn)
    #DeltaEtot_t_i    = np.array(DeltaEtot_t_i)
    #print "Dumping to .pkl file..."
    #cPickle.dump( (DeltaEmodes_i_tn, DeltaEtot_t_i), open(modes_pkl,"w") );
    #
    #print "TrackModes.py complete."

if __name__ == "__main__":
	main()
