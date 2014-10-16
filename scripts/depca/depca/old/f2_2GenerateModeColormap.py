import numpy as np
import matplotlib.pylab as plt
import h5py
from f1SidechainCorr import h5crtag
import f2_1TrackModes as TrackModes
import ConfigParser
import argparse
import colorbrewer
from scipy.interpolate import interp1d

import math

def FMO_site2sidechain(FMOgro):
    count = 7
    count_BCX = 0
    FMO_sidechain_map = {}
    FMO_chromo_map    = {}
    for l in FMOgro.split('\n'):
        larr = l.split()
        if len(larr) == 9 and len(larr[0]) <= 6:
            i = int(larr[0][:-3])
            res = larr[0][-3:]
            if not i in FMO_sidechain_map.values():
                print res,
                if res == "BCL":
                    FMO_sidechain_map[count_BCX] = i
                    FMO_chromo_map[count_BCX + 1] = i
                    count_BCX += 1
                else:
                    FMO_sidechain_map[count] = i
                    count += 1
    print FMO_sidechain_map
    print FMO_chromo_map
    return FMO_sidechain_map, FMO_chromo_map

def rgb2hexcode(rgb):
    assert(len(rgb)==3)
    hexcode = "#"
    for c in rgb:
        twodigit = hex(int(c))[2:]
        if len(twodigit) > 2:
            print twodigit
        if len(twodigit) == 1:
            twodigit = "0"+twodigit
        hexcode+=twodigit.upper()
    return hexcode

def InterpCode(x, diverging_scale):
    xscale = np.linspace(-1,1,len(diverging_scale))
    diverging_scale = np.array(diverging_scale)
    f = interp1d(xscale, diverging_scale, axis=0, kind = 'cubic')
    return rgb2hexcode(f(x))



def GenerateModeColormap_ChimeraPrint(mode, site_num, numscale, FMO_sidechain_map, FMO_chromo_map):
        BrewerScale = colorbrewer.PiYG
        BrewerScale = [rgb2hexcode(a) for a in BrewerScale[numscale]]
	half_scale = float((numscale - 1) / 2)
        out = ""
       
        assert numscale % 2 == 1
        assert numscale <= 11

	out += "# Site Number is: {}\n".format(site_num + 1)
	out += "# Max Value: {}\n".format(max(mode))
	out += "# Min Value: {}\n".format(min(mode))

	max_val = max(abs(max(mode)), abs(min(mode)))
        print mode.shape
	for i, x in enumerate(mode):
		#sigmoid_val = (x / max_val) * 20.
		#new_val = math.erf(sigmoid_val) / math.erf(2)
		#val =   x / max_val * half_scale
                #if val >= .5 or val <= -.5:
                #    print val, max_val
                #color_num = int(round(val))
		#out += "#{}\n".format(x/max_val)
		#assert color_num >= -half_scale and color_num <=half_scale
		#color_code = BrewerScale[color_num + int(half_scale)]
                color_code = InterpCode(x / max_val, colorbrewer.PiYG[numscale])

		sidechain= FMO_sidechain_map[i]

		out += "color {} :{}\n".format(color_code,sidechain)

	for site_i in xrange(7):
		chromo = FMO_chromo_map[site_i + 1]

		if site_i == site_num:
			out += "color #FF00FF :{}\n".format(chromo)
		else:
			out += "color {} :{}\n".format(BrewerScale[int(half_scale)], chromo)

	out += "# Impact, (i.e. sum) is: {}".format(sum(mode))
        return out



def main():
    config = ConfigParser.RawConfigParser()
    config.read('./f0postProcess.cfg')
    h5file    = config.get('sidechain','h5file')
    corrtag   = config.get('sidechain','corr_h5tag')
    num_scale = config.getint('chimera','num_scale')
    gro_file  = config.get('chimera','gro_file')

    parser = argparse.ArgumentParser(description="Generate a colormap for Chimera that visualizes the modes that come from the PCA matrix stored in the HDF file")
    parser.add_argument("site", type=int, help="Site that the data is requested for, use 1-based indexing. No error checking.")
    parser.add_argument("mode", type=int, help="Modes to generate outputs for")
    parser.add_argument("dest", type=str, default="modes.com", help="location to save output.")
    parser.add_argument("-debug", action="store_true")
    args = parser.parse_args()

    site_num = args.site - 1
    mode_num = -args.mode

    with open(gro_file,'r') as f:
        FMOgro = f.read()
    FMO_sidechain_map, FMO_chromo_map = FMO_site2sidechain(FMOgro)

    #help(colorbrewer)
    
    if args.debug:
        BrewerScale = colorbrewer.PiYG[11]
        x = .1
        InterpCode(x, BrewerScale)
        x = .2
        InterpCode(x, BrewerScale)
        exit(0)

    

    with h5py.File(h5file,'r') as f:
        corr   = f[corrtag+h5crtag]
	vi, wi, impact_i = TrackModes.ComputeModes(corr)
	out = GenerateModeColormap_ChimeraPrint(wi[site_num][mode_num], site_num, num_scale, FMO_sidechain_map, FMO_chromo_map)
        site_amp = sorted(wi[site_num][mode_num], key=abs, reverse=True)
        print site_amp[0:4]
    with open(args.dest, 'w') as f:
        f.write(out)

if __name__ == "__main__":
	main()
