import argparse
import h5py
import numpy as np
import re


# Specification: Take (input1) [csv format] and append into (input2) [hdf5 format] under database name (input3)
# Note: Should be chained into a bash script to handle multiple files
# Note: Should have a delete (y/n) option

def stripnum(fname):
    pre_name = fname.split('.')[-2]
    y = re.search('\d*$', pre_name)
    if y:
        print pre_name[y.start():]
        return int(pre_name[y.start():])
    else:
        raise RuntimeError("Cannot find numeric value for search")


def main():
    print "~~~csv2hdf5~~~~"
    parser = argparse.ArgumentParser(description='Take (input1-N) [csv format] and append into (input2) [hdf5 format] under database name (input3). Requires that .csv file fits in memory. Also, requires that csv row width is 357 (value from FMO).')
    parser.add_argument('hdf_file', type=str, help=".hdf5 destination file.")
    parser.add_argument('hdf_dsname', type=str, help="database entry name within hdf_file")
    parser.add_argument('csv_file', type=str, nargs='*', help=".csv file(s) to merge into hdf5.")
    parser.add_argument('--open_with_write', action='store_true', help="when set, delete or create database at hdf_file ('w' flag); without flag, append or create database at hdf_file ('a' flag).")
    parser.add_argument('--sort', action='store_true', help="when set, sort .csv files by their numbers just preceeding their final extension; otherwise, use the order passed into the function.")
    parser.add_argument('-frames_per_file', type=int, default=50000, help="number of frames in each file; must be correct for script to function. (Default = 50,000)")
    parser.add_argument('-dt', type=float, default=None, help="timestep, in picoseconds")

    args = parser.parse_args()
    open_flag = "a"
    if args.open_with_write:
        open_flag = "w"

    if args.sort:
        args.csv_file.sort(key=stripnum)


    print "Writing to {}/{} with '{}'".format(args.hdf_file, args.hdf_dsname, open_flag)

    with h5py.File(args.hdf_file,open_flag) as h5_out:
        dsshape = (None, 7, 357)
        # Check if the dataset exists
        h5keys = h5_out.items()
        goodkeys = [key[0] == args.hdf_dsname for key in h5keys]
        if any(goodkeys):
            ds = h5_out[args.hdf_dsname]
        else:
            ds = h5_out.create_dataset(args.hdf_dsname, shape=(0,7,357), maxshape=dsshape)

        if args.dt:
            ds.attrs['dt'] = args.dt

        for csv_file in args.csv_file:
            print "Loading .csv file into memory: {}, {} rows....".format(csv_file, args.frames_per_file)
            with open(csv_file, 'r') as csv_in:
                f = csv_in.read().strip()
                #csv_dat.append(np.array(l.split()))
            l_arr = f.split("\n")
            l_arr = np.array([[float(x) for x in l.split(',')] for l in l_arr])
            x = l_arr[1,:]
            print "Reshaping csv data to ({},{},{})...".format(args.frames_per_file, 7,357)
            l_arr.shape = (args.frames_per_file, 7, 357)
            y = l_arr[0,1,:]
            assert all(y == x)

            oldlen = ds.shape[0]
            newlen = oldlen + l_arr.shape[0]
            
            print "Resizing hdf5 dataset: ({},{},{}) -> ({},{},{})".format(oldlen, dsshape[1],dsshape[2], newlen, dsshape[1], dsshape[2])
            ds.resize(newlen, axis=0)
            ds[oldlen:newlen,:,:] = l_arr[:]



if __name__ == "__main__":
    main()
