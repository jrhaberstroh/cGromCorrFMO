import numpy as np
import cPickle
import scipy.linalg as LA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import math
import sys
from copy import deepcopy
from SidechainCorr import SidechainRead
import ConfigParser



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


def PlotLogSpectrum(vi, impact_i, floor = 1E-15):
	Ni = len(vi)
	for i in xrange(Ni):
		vlog = [math.log( abs((x * imp)) + floor, 10) for x, imp in zip(vi[i], impact_i[i])]
		plt.plot(vlog, 'ro')
		plt.plot((math.log(floor,10),) * len(vi[i]))
		plt.savefig('plot_spectrum.pdf')
		plt.savefig('plot_spectrum.png')
		#plt.show()



def DeltaModeTracker(E_t_ij, modes_inj):
	"""
	DeltaModeTracker: Takes E(t) np.ndarray and an iterable of modes, and plots the deviation of dE(t) [after taking the average for each site]
	                  as contributed from the modes, as compared to the actual dE(t)

	\input 	E_t_ij 		- a collection of Eij as a function of time
		modes_inj 	- a collection of n normalized modes (in j coordinates) for each site i

	\output DeltaEm_t_in	- a collection of the variations in energetics for the n modes for site i as a function of time
	"""

	print "Running DeltaModeTracker..."
	E_t_ij = np.array(E_t_ij)

	DeltaEmodes_t_in = []
	DeltaEtot_t_i = []

	Ntimes = len(E_t_ij)
	Nsites = len(E_t_ij[0])
	Ncoarse = len(E_t_ij[0][0])

	E_avg_ij = np.zeros( (Nsites, Ncoarse) )
	print "Ntimes: ", Ntimes
	print "Nsites: ", Nsites
	print "Ncoarse: ", Ncoarse
	print "Averages array size: ", E_avg_ij.shape
	
	print "Computing time-average for each site..."
	for i in xrange(Nsites):
		for j in xrange(Ncoarse):
			x = E_t_ij[:,i,j]
			E_avg_ij[i][j] = np.mean(E_t_ij[:,i,j])


#	print "Compute dE on each mode, at each time..."
#	for t, Eij in enumerate(E_t_ij):
#		if t%100 == 0:
#			print "t =",t
#		DeltaEmodes_t_in.append([])
#		DeltaEtot_t_i.append([])
#		for i, Ej in enumerate(Eij):
#			DeltaEmodes_t_in[t].append([])
#			for n, Mj in enumerate(modes_inj[i]):
#				weight_n = sum(Mj)
#				DeltaEmodes_t_in[t][i].append( (np.dot(Ej, Mj) - np.dot(E_avg_ij[i], Mj)) * weight_n )
#			DeltaEtot_t_i[t].append(sum(Ej))

	# First, check that the modes are normalized, then weight the modes
	print "Checking normalization and applying weights..."
	for m_nj in modes_inj:
		for m_j in m_nj:
			assert abs(LA.norm(m_j) - 1) < 1E-6
			m_j *= sum(m_j)
	
	#Then, vectorize the code:
	print "Running vectorized mode tracker computation..."
	dEmode_i_tm = []
	for i in xrange(Nsites):
		y = np.matrix(np.inner(E_avg_ij[i,:], modes_inj[i,:,:]))
			# y is (1,m)
		x = np.matrix(np.inner(E_t_ij[:,i,:], modes_inj[i,:,:]))
			# x is (t,m)

		# Matrix subtraction will repeat rows of y to subtract from x
		dEmode_i_tm.append(np.array(x - y))


	print "Running total dE(t) computation..."
	#   [1 x i] matrix
	y = np.matrix(np.inner(E_avg_ij, np.ones(Ncoarse)))
	#   [t x i] matrix           
	x = np.matrix(np.inner(E_t_ij, np.ones(Ncoarse)))
	
	dEtot_t_i = np.array(x - y)

	print "Done."

	return dEmode_i_tm, dEtot_t_i








#--------------------------------------------------------------


def main():
	config = ConfigParser.RawConfigParser()
	config.read('./.postProcess.cfg')
	filename = config.get('sidechain','pkl_file')
	csv_filename = config.get('sidechain','csv_file')
	
	modes_pkl = config.get('modes','modes_pkl')
	resave_modes = config.getboolean('modes','resave_modes')
	plot_spectrum = config.getboolean('modes','plot_spectrum')
	time_max = config.getint('modes','time_max')

	
	print "Loading "+ filename + "..."
	corr, Avg_Eij = cPickle.load(open(filename, 'r'))
	print "\tLoaded"
	

	print "Computing Modes..."
	vi, wi, impact_i = ComputeModes(corr)
	print "\tModes Computed!"

	if plot_spectrum:
		PlotLogSpectrum(vi, impact_i)

	if not resave_modes:
		print "resave_modes == no, exiting TrackModes.main()"
		exit(101)
	else:
		print "resave_modes = yes, saving the computed trajectories for PlotModeTimeseries..."

	N = time_max
	E_t_i = []
	nread = 0

	print "Reading",N,"times of sidechain data from "+csv_filename+" to plot timeseries..."
	E_t_i, nread = SidechainRead(csv_filename, N)
	print "\tActually read", nread, "times from the file."

	DeltaEmodes_i_tn, DeltaEtot_t_i = np.array(DeltaModeTracker(E_t_i, wi))
	DeltaEmodes_i_tn = np.array(DeltaEmodes_i_tn)
	DeltaEtot_t_i    = np.array(DeltaEtot_t_i)
	print "Dumping to .pkl file..."
	cPickle.dump( (DeltaEmodes_i_tn, DeltaEtot_t_i), open(modes_pkl,"w") );
	
	print "TrackModes.py complete."

if __name__ == "__main__":
	main()
