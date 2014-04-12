import csv
import cPickle
import matplotlib.pylab as plt
import numpy as np
from copy import deepcopy


def SidechainPlot(filename, cutoff = 3):
	dE_t_i = []
	counter = 0;
	times = 0;
	with open(filename,'r') as f:
		read = csv.reader(f)

		dE_i = []
		myArr = []

		for pos, csvline in enumerate(read):
			myArrOld = deepcopy(myArr)
			myArr = [float(x) for x in csvline]
			dE_i.append(myArr)
			if times >= cutoff:
				break
			if times < cutoff:
				print times,", ",counter
				if (len(myArrOld) > 0):
					plt.plot(np.array(myArrOld) - np.array(myArr))
					plt.show()
				else:
					plt.plot(myArr)
					plt.show()
			counter += 1
			if (counter == 7):
				counter = 0;
				times = times + 1;
				dE_t_i.append(np.array(deepcopy(dE_i)))
				dE_i = []



if __name__ == "__main__":
	SidechainPlot("./fmoGap_t3k_all7.csv")
