import argparse
import h5py
import numpy as np


# Specification: Take (input1) [csv format] and append into (input2) [hdf5 format] under database name (input3)
# Note: Should be chained into a bash script to handle multiple files
# Note: Should have a delete (y/n) option

def main():
    print "---------------------------------csv2hdf5, hurrah!---------------------------------------"
    parser = argparse.ArgumentParser(description='Take (input1-N) [csv format] and append into (input2) [hdf5 format] under database name (input3). Requires that .csv file fits in memory.')
    parser.add_argument('csv_file', type=str, nargs='+', help=".csv file(s) to merge into hdf5")
    parser.add_argument('hdf_file', type=str, help=".hdf5 destination file.")
    parser.add_argument('hdf_dsname', type=str, help="database entry name within hdf_file")
    parser.add_argument('--nodelete', action='store_true', help="when set, append or create database at hdf_file; otherwise, delete or create database at hdf_file.")

    args = parser.parse_args()
    open_flag = "a"
    if not args.nodelete:
        open_flag = "w"
    else:
        raise NotImplementedError("Program does not know how to run without --nodelete; maybe you could write this patch?")


    with h5py.File(args.hdf_file,open_flag) as h5_out:
        dsshape = (None, 7, 357)
        if not args.nodelete:
            ds = h5_out.create_dataset(args.hdf_dsname, shape=(0,7,357), maxshape=dsshape, chunks=(400, 7, 357))

        for csv_file in args.csv_file:
            with open(csv_file, 'r') as csv_in:
                f = csv_in.read().strip()
                #csv_dat.append(np.array(l.split()))
            l_arr = f.split("\n")
            l_arr = np.array([[float(x) for x in l.split(',')] for l in l_arr])
            x = l_arr[1,:]
            l_arr.shape = (400, 7, 357)
            y = l_arr[0,1,:]
            assert all(y == x)

            oldlen = ds.shape[0]
            newlen = oldlen + l_arr.shape[0]
            ds.resize(newlen, axis=0)
            print ds.shape
            ds[oldlen:newlen,:,:] = l_arr[:]
            
        





if __name__ == "__main__":
    main()
