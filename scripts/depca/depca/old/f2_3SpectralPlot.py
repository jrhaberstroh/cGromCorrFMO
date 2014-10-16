import ConfigParser
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anime
import pandas as pd
import argparse

def CtAnimate(Ct_ab):
    fig = plt.figure()
    ims = []
    #im = plt.imshow(np.log(Ct_ab[0,:,:]))
    for C_ab in Ct_ab[0:100,:,:]:
        im = plt.imshow(np.log(C_ab))
        im.set_cmap('gray')
        label = plt.text(20,20, "t = {}".format(np.amax(C_ab)),  bbox=dict(facecolor='red', alpha=0.5))
        ims.append([im, label])
    ani = anime.ArtistAnimation(fig,ims, interval=20,repeat_delay=0)
    plt.show()


def main():
    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')
    h5_filename = config.get('sidechain','h5file')
    h5time = config.get('sidechain','time_h5tag')
    h5corr = config.get('sidechain','corr_h5tag')
    h5ct = config.get('sidechain','ct_h5tag')

    parser = argparse.ArgumentParser(description = "Module to visualize temporal correlations")
    parser.add_argument('-chromo', type=int, default=1, help="Chromophore to compute correlations for")
    args = parser.parse_args()

    h5ct = "".join([h5ct, "{}".format(args.chromo)])

    with h5py.File(h5_filename,'r') as f:
        Ct_ab = f[h5ct]
        CtAnimate(Ct_ab)
    


if __name__=="__main__":
    main()
