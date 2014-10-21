import depca.visualize.pcaviz as pca
import depca.dedata as dedata




def PlotGapDist():
    with dedata.dEData('./f0postProcess.cfg') as data:
        print "Running with {}...".format(type(data))
        sidechain = data.GetModes_hdf()
        print "Plotting..."
        legend = None
        for i in xrange(1,10):
            subsamples = pca.GetSubsamples(sidechain[:,0,i], .05, 10)
            print subsamples.shape
            print subsamples[0,:]
            pca.Plot1DHist(subsamples, legend=legend, displace_by=0.0)


if __name__=="__main__":
    PlotGapDist()
