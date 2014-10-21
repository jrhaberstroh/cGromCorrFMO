import numpy as np
import scipy.linalg as LA
import matplotlib.pylab as plt
import h5py
from copy import deepcopy
import ConfigParser
import argparse
import random

def GetSubsamples(array1d, fraction, quantity):
    num_in_sample = int(len(array1d) * fraction)
    subsamples = np.zeros((quantity, num_in_sample))
    print len(array1d)
    print num_in_sample
    for i in xrange(quantity):
        sample_ind = random.sample(xrange(len(array1d)), num_in_sample)
        print array1d[sample_ind]
        subsamples[i,:] = array1d[sample_ind]
    return subsamples

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

def Plot1DHist(entries, residual = None, displace_by = 1.0, plottype = 'display', fname = "plot2d", legend= [], free_energy=True, parabola=False):
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

    if len(entries.shape) == 1:
        Plot1DHistMode(entries, offset=0, plots=plots, free_energy=free_energy, parabola=False)
    
    else:
        for i in xrange(entries.shape[0]):
            Plot1DHistMode(entries[i,:], offset=i*displace_by, plots=plots, free_energy=free_energy, parabola=False)
    
    if residual:
        Plot1DHistMode(residual, offset=-1, plots=plots, free_energy=free_energy, parabola=False)

    print plots
    print legend
    if legend and len(legend) > 0:
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



def ComputeCt_FFT(f1_t, f2_t, corr_len=None):
    if len(f1_t) != len(f2_t):
        raise TypeError("ComputeCt trying to correlate functions of different sizes: {} and {}".format(len(f1_t), len(f2_t)))
    if not corr_len:
        corr_len = len(f1_t)
    numavg  = len(f1_t) - np.arange(corr_len)
    mycorr = np.real(np.fft.ifft( np.fft.fft(f2_t) * np.conj(np.fft.fft(f1_t)) )[0:corr_len])/numavg
    return mycorr




