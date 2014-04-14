import numpy as np
import cPickle
import scipy.linalg as LA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import math
import sys
from copy import deepcopy
from f1SidechainCorr import SidechainRead
import ConfigParser

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



def PlotScatterFirstTwo(dEmodes_i_tn, dEtot_t_i, mode1, mode2, i = 0, plotType = 1):
	assert(type(dEmodes_i_tn) == np.ndarray)
	assert(mode1 <= dEmodes_i_tn.shape[2])
	assert(mode2 <= dEmodes_i_tn.shape[2])
	Ntimes = dEmodes_i_tn.shape[1]

	mode1 = -mode1 - 1
	mode2 = -mode2 - 1
	x = dEmodes_i_tn[i,:,mode1]
	y = dEmodes_i_tn[i,:,mode2]
	
	if (plotType == 1):
		plt.scatter(x, y)
		plt.title('Modes ['+str(-mode1)+"] and ["+str(-mode2)+"] on site "+str(i+1)+" (all numbers 1-based indexing)")
		print "Saving scatterplot figure..."
		plt.savefig("plot_sctr_"+str(-mode1)+"v"+str(-mode2)+"s"+str(i+1)+".png")
		plt.savefig("plot_sctr_"+str(-mode1)+"v"+str(-mode2)+"s"+str(i+1)+".pdf")
		#plt.show()

	if (plotType == 2):
		xwidth = (max(x) - min(x)) / 2. * 1.10
		xmid   = (max(x) + min(x)) / 2.
		ywidth = (max(y) - min(y)) / 2. * 1.10
		ymid   = (max(y) + min(y)) / 2.
		xedges = np.linspace(xmid - xwidth, xmid + xwidth, 50)
		yedges = np.linspace(ymid - ywidth, ymid + ywidth, 50)
	
		print "Saving histogram figure..."
		H, xedges, yedges = np.histogram2d(x, y, bins=(xedges,yedges))
		im = plt.imshow(H, interpolation='nearest', origin='low', extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
		plt.title('Modes ['+str(-mode1)+"] and ["+str(-mode2)+"] on site "+str(i+1)+" (all numbers 1-based indexing)")
		plt.savefig("plot_hist_"+str(-mode1)+"v"+str(-mode2)+"s"+str(i+1)+".png")
		plt.savefig("plot_hist_"+str(-mode1)+"v"+str(-mode2)+"s"+str(i+1)+".pdf")
		#plt.show()


def main():
	config = ConfigParser.RawConfigParser()
	config.read('./.postProcess.cfg')
	filename = config.get('sidechain','pkl_file')
	
	modes_pkl = config.get('modes','modes_pkl')

	plot_mode = config.get('plotter','plot_mode')
	plot_type = config.getint('plotter','plot_type')
	num_modes = config.getint('plotter','num_modes')
	site_num = config.getint('plotter','site_num') - 1
	mode1 = config.getint('plotter','scatter_1') - 1
	mode2 = config.getint('plotter','scatter_2') - 1

	save_plot = config.getboolean('plotter','save_to_disc')

	print "plot_modes_only is on!"
	print "Loading modes from modes_pkl =", modes_pkl+"..."
	dEmodes_i_tn, dEtot_t_i = cPickle.load(open(modes_pkl, 'r'))

	dEmodes_i_tn = np.array(dEmodes_i_tn);
	dEtot_t_i = np.array(dEtot_t_i);

	print "Plotting modes from modes_pkl =", modes_pkl+"..."
	print "Array shape, modes: ", dEmodes_i_tn.shape
	print "Array shape, tot: ", dEtot_t_i.shape
	print "Num times: ", dEmodes_i_tn.shape[1]


	if plot_mode == "separate":
		dEModesSeparate(dEmodes_i_tn, dEtot_t_i, num_modes, site_num)
	if plot_mode == "total":
		PlotDeltaTimeseries(dEmodes_i_tn, dEtot_t_i, num_modes, site_num)
	if plot_mode == "scatter":
		print "Beginning scatterplot figure..."
		PlotScatterFirstTwo(dEmodes_i_tn, dEtot_t_i, mode1, mode2, site_num, plot_type)

	print "\tDone!"


	sys.exit(1)




if __name__ == "__main__":
	main()
