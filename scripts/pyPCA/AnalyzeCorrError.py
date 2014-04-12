import numpy as np
import ConfigParser
import h5py
import SidechainCorr
import matplotlib.pylab as plt
import scipy.linalg as la

def CorrCompare(Corr1, Corr2):
	return 0

def PlotMtxError(Corr_w):
	max_val = 1
	min_val = -.1

	AvCorr = np.sum(Corr_w, axis=0)
	dCorr = Corr_w - AvCorr
	errCorr = np.log10(np.sqrt(np.einsum("i...,i...", dCorr, dCorr)) / np.absolute(AvCorr) / np.sqrt(Corr_w.shape[0]))
	#print errCorr.shape
	#print errCorr

	plt.rcParams.update({'font.size':6, 'font.weight':'bold'})
	for i in xrange(errCorr.shape[0]):
		plt.subplot(2,7,i + 1)
		plt.title("SITE "+str(i+1)+":: \nHistogram of errors in corr. mtx.")
		plt.hist(errCorr[0,:,:].flatten(), 256, range=(min_val, max_val))
		plt.xlabel("log_10(sigma)")
		plt.ylabel("Count")
	
		plt.subplot(2,7,i + 7 +1)
		plt.imshow(errCorr[0,:,:], vmin = min_val, vmax = max_val)
		cbar = plt.colorbar(shrink=.25, aspect=40)
		cbar.set_label("log_10(sigma)")
		plt.set_cmap('gist_yarg')
		plt.title("SITE "+str(i+1)+":: \nError in corr. matx. values")
		plt.xlabel("Site i")
		plt.ylabel("Site j")
	plt.show()


def PlotEvDrift(Corr_tot, Corr_w):
	Nsite = Corr_tot.shape[0]
	v_w = [la.eigh(Corr_tot[i]) for i in xrange(Nsite)]
	V = [v_w[i][0][::-1] for i in xrange(Nsite)]
	W = [v_w[i][1][::-1] for i in xrange(Nsite)]
	plt.plot(np.log10(V[0]), 'ro')
	plt.show()



if __name__ == "__main__":
	config = ConfigParser.RawConfigParser()
	config.read('./.postProcess.cfg')
	csv_filename = config.get('sidechain','csv_file')
	pkl_filename = config.get('sidechain','pkl_file')
	h5_filename = config.get('sidechain','h5_file')
	h5_tag = config.get('sidechain','h5_tag')

	with h5py.File(h5_filename,'r') as f:
		ds_name = h5_tag+SidechainCorr.h5detag
		print "Loading "+ds_name+" from file "+h5_filename+"..."
		dedset = f[ds_name];
		print "\tShape of data", dedset.shape

		print "\tPlotting Example data..."
		for i in xrange(7):
			plt.subplot(3,3,i+1)
			plt.plot(np.sum(dedset[:,i,:], axis=1) * 350)
			plt.title("SITE "+str(i+1))
			plt.ylabel("Energy, cm-1")
		plt.show()

		num_times = tsdset.shape[0]
		
		partial = 5;
		windows = 10;

		step_size   = num_times / windows;
		window_size = num_times / partial;
		
		start = [ ((num_times - window_size) / windows) * i for i in np.arange(windows+1) ]
		end   = [ st + window_size for st in start ]
		
		print "Windows to use: ", zip(start,end)

		types = ['Correlation_iab','Average_ia']

		Corr_w = []
		Avg_w = []

		Corr_tot = SidechainCorr.AvgAndCorrelateSidechains(\
				tsdset, tsdset.shape[0])[0]

		for (st,ed) in zip(start,end):
			Corr_new, Avg_new = SidechainCorr.AvgAndCorrelateSidechains(\
					tsdset[st:ed,:,:], window_size)
			Corr_w.append(Corr_new)
			Avg_w.append(Avg_new)

		Corr_w = np.array(Corr_w)
		Avg_w = np.array(Avg_w)
	
		print Corr_w.shape
		print Avg_w.shape
		
		#print "Estimating error on each correlation element..."
		#PlotMtxError(Corr_tot, Corr_w)

		print "Plotting drift in eigenvalues relative to best eigenvalue"
		PlotEvDrift(Corr_tot, Corr_w)
		

		#print CorrAndAvg_w_type[0][1].shape


