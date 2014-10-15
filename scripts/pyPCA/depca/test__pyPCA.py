import f1SidechainCorr as SidechainCorr
import f2_1TrackModes as TrackModes
import f3PlotModeTimeseries as PlotModeTimeseries
import numpy as np
import random



n_times = 1000
random.seed(90210)

E_t_ij = np.array((np.array((np.array([0.,0.,0.,0.]),) * 1),) * n_times)

for t in xrange(n_times):
	r = random.gauss(0,1)
	#print r
	E_t_ij[t,0,1] = 1. + r
	E_t_ij[t,0,3] = 1. + r


#print E_t_ij

Corr, Avg = SidechainCorr.AvgAndCorrelateSidechains(E_t_ij, 10000)

print Corr
print Avg
print type(Corr)
print type(Avg)

vi, wi, impact_i = TrackModes.ComputeModes(Corr)

print vi
print wi
print impact_i

xi = []
for i in xrange(len(vi)):
	xi.append([x * y for x, y in zip(vi[i],impact_i[i])])

print "TOTAL IMPACTS:"
print xi

#TrackModes.PlotLogSpectrum(vi, impact_i)

dEmodes_i_tm, dEtot_t_i = TrackModes.DeltaModeTracker(E_t_ij, wi)
dEmodes_i_tm = np.array(dEmodes_i_tm);
dEtot_t_i = np.array(dEtot_t_i);

print dEmodes_i_tm.shape

PlotModeTimeseries.PlotDeltaTimeseries(dEmodes_i_tm, dEtot_t_i, 2)
