import numpy as np
import csv
import cPickle
import h5py
import sys
import matplotlib.pylab as plt
import ConfigParser
from copy import deepcopy

h5tstag = "timeseries_tia"
h5crtag = "corr_ij"
h5eavtag = "eav_ij"
h5detag = "deltaseries_tia"
h5mstag = "modeseries_tia"


def SidechainRead(filename, cutoff = 3):
	dE_t_i = []
	counter = 0;
	times = 0;
	with open(filename,'r') as f:
		read = csv.reader(f)

		dE_i = []
		myArr = []

		for pos, csvline in enumerate(read):
			myArrOld = myArr
			myArr = np.array([float(x) for x in csvline])
			dE_i.append(myArr)
			if times >= cutoff:
				break
		#	if times < cutoff:
		#		print times,", ",counter
		#		if (len(myArrOld) > 0):
		#			plt.plot(np.array(myArrOld) - np.array(myArr))
		#			plt.show()
		#		else:
		#			plt.plot(myArr)
		#			plt.show()
			counter += 1
			if (counter == 7):
				counter = 0;
				times = times + 1;
				dE_t_i.append(np.array(dE_i))
				dE_i = []
	return np.array(dE_t_i), times

def SidechainReadv2(filename, t_start = 0, t_end = 3):
	dE_t_i = []
	counter = 0;
	times = 0;
	with open(filename,'r') as f:
		read = csv.reader(f)

		dE_i = []
		myArr = []

		for pos, csvline in enumerate(read):
			if times >= t_start:
				myArrOld = myArr
				myArr = np.array([float(x) for x in csvline])
				dE_i.append(myArr)
			if times >= t_end:
				break
		#	if times < cutoff:
		#		print times,", ",counter
		#		if (len(myArrOld) > 0):
		#			plt.plot(np.array(myArrOld) - np.array(myArr))
		#			plt.show()
		#		else:
		#			plt.plot(myArr)
		#			plt.show()
			counter += 1
			if (counter == 7):
				if times >= t_start:
					dE_t_i.append(np.array(dE_i))
					dE_i = []
				counter = 0;
				times = times + 1;
	return np.array(dE_t_i), times

def AvgAndCorrelateSidechains(E_t_ia, endTime = 3):
	endTime = min(endTime, len(E_t_ia))

	AvgEia    = E_t_ia[0:endTime,:,:].sum(axis=0)
	AvgEiaEib = np.array( [ np.tensordot(E_t_ia[0:endTime,i,:], E_t_ia[0:endTime,i,:], axes=(0,0)) for i in xrange(E_t_ia.shape[1]) ] )
	
	AvgEia    /= endTime
	AvgEiaEib /= endTime

	Corr_i_ab = AvgEiaEib - np.einsum('...j,...k',AvgEia,AvgEia)

	return Corr_i_ab, AvgEia


def ShowData(E_t_i):
	cutoff = 2

	for t, Et in enumerate(E_t_i):
		print t
		if t >= cutoff:
			break
		for i, Ei in enumerate(Et):
			print t,", ", i
			plt.plot(Ei)
			plt.show()

def main():
	config = ConfigParser.RawConfigParser()
	config.read('./.postProcess.cfg')
	csv_filename = config.get('sidechain','csv_file')
	pkl_filename = config.get('sidechain','pkl_file')
	h5_filename = config.get('sidechain','h5_file')
	h5_tag = config.get('sidechain','h5_tag')
	t_start = config.getint('sidechain','t_start')
	t_end = config.getint('sidechain','t_end')

	times = 10000
	if (t_end < 0):
		# map t_end cyclically around "times", which is really not the end of the array...
		t_end = ((((times + t_end + 1) % times) + times) % times)

	print "Reading sidechains from "+csv_filename+"..."
	E_t_ia, times_read = SidechainReadv2(csv_filename, t_start, t_end)
	E_t_ia = np.array(E_t_ia)
	print "\tRead",len(E_t_ia),"lines."
	
	print "Computing Correlations..."
	corr_iab,Avg_Eia = AvgAndCorrelateSidechains(E_t_ia, times)
	#ShowData(x)

	openflag = 'a'
	opendesc = {'r':"Readonly, file must exist", 'r+':"Read/write, file must exist", \
			'w':"Create file, truncate if exists", 'w-':"Create file, fail if exists", \
			'a':"Read/write if exists, create otherwise"}
	print "Opening database "+h5_filename+" with option "+openflag+" [i.e. "+opendesc[openflag]+"]..."
	with h5py.File(h5_filename,openflag) as f:
		try:
			print "\tSaving to "+h5_tag+h5tstag+" timeseries of shape",E_t_ia.shape,"..."
			tsdset = f.require_dataset(h5_tag+h5tstag, E_t_ia.shape, dtype='f');
			tsdset[...] = E_t_ia
			print "\tSaving to "+h5_tag+h5crtag+" correlation matrix of shape",corr_iab.shape,"..."
			crdset = f.require_dataset(h5_tag+h5crtag, corr_iab.shape, dtype='f');
			crdset[...] = corr_iab
			print "\tSaving to "+h5_tag+h5eavtag+" averages of shape",Avg_Eia.shape,"..."
			eavdset = f.require_dataset(h5_tag+h5eavtag, Avg_Eia.shape, dtype='f');
			eavdset[...] = Avg_Eia
	
			# Create the delta-timeseries
			dE_tia = E_t_ia - Avg_Eia
			print "\tSaving to "+h5_tag+h5detag+" timeseries of shape",dE_tia.shape,"..."
			dedset = f.require_dataset(h5_tag+h5detag, dE_tia.shape, dtype='f');
			dedset[...] = dE_tia
		except TypeError:
			print "\n\nERROR: Error in write encountered, you should consider changing the h5_tag to a value that has not been used."
			raise

	print "Database closed"

	#cPickle.dump((corr,Avg_Eij), open(pkl_filename,"wr"))


if __name__ == "__main__":
	main()
