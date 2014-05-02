import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import argparse
import h5py

def main():
    parser = argparse.ArgumentParser(description = "Generate SHO spectral density bath data to be used as proxy data for understanding the results of pyPCA module analysis tools.")
    parser.add_argument('-nSHO', type=int, default=357, help="Number of harmonic oscillators to use in the dynamics")
    parser.add_argument('-dt', type=float, default=.01, help="Timestep for simulation, ps. Stored in hdf5 file.")
    parser.add_argument('-T', type=float, default=2000., help="Total time to generate, ps.")
    parser.add_argument('-h5out', default='SHO.hdf5', help="Location for the output of the dynamics")
    parser.add_argument('-h5dset', default='alltimes', help="Name of the timeseries in the hdf5 database")
    parser.add_argument('-DOF', type=int, default=1, help="Number of DOF merged into each time slice")
    args = parser.parse_args()
    
    n_steps = round(args.T / args.dt) + 1
    nSHO = args.nSHO
    rand.seed(90211)
    #amp = np.array([np.exp(np.log(8) * abs(rand.normal())) for x in [0]*357])
    omega_cm = np.linspace(0,1500,5000)
    beta = 1
    freq = 1 * np.tanh(beta * omega_cm / 2) *  ( ( 100 / (20000 + np.square(omega_cm)))  + \
                                                 ( 100 / (30000 + np.square(omega_cm)))  + \
                                                 ( 1000 / (2 * (10000 + np.square(omega_cm - 1000)))) + \
                                                 ( 1000 / (2 * (10000 + np.square(omega_cm - 500 ))))  )
    freq /= sum(freq)
    CDF_W = np.cumsum(freq)

    planck_cmps = 5.308
    A  = np.linspace(0.,100., 101)
    pA = np.reciprocal(A)
    pA[0]  = pA[1]
    CDF_A  = np.cumsum(pA / np.sum(pA))
    #aSHO     = rand.gamma(10., .1, (args.DOF, nSHO))
    aSHO     = np.interp(rand.rand(args.DOF, nSHO), CDF_A, A)
    pSHO_rad = rand.rand(args.DOF,nSHO) *2.*np.pi
    wSHO_ps  = np.interp(rand.rand(args.DOF, nSHO), CDF_W, omega_cm/planck_cmps) 

    H, bins = np.histogram(wSHO_ps, weights=aSHO, bins=100)
    binsctr = .5 * bins[1:] + .5 * bins[:-1]

    plt.plot(binsctr, H/float(np.sum(H * (bins[1:] -  bins[:-1]))))
    plt.plot(omega_cm/planck_cmps, freq / float(omega_cm[1] - omega_cm[0])*planck_cmps)
    plt.xlabel("1/ps, frequency")
    plt.ylabel("Density of SHO")
    plt.title("Spectral density vs actual modes")

    plt.show()

    t = np.linspace(0,args.T, n_steps)

    if False:
        with h5py.File(args.h5out,'w') as f:
            E_tij = f.create_dataset(args.h5dset, (n_steps, 1, nSHO), fillvalue=0.)
            E_tij.attrs['dt'] = args.dt
            for i in xrange(args.DOF):
                E_tij[:] += np.cos(np.outer(t,wSHO_ps[i,:]) + pSHO_rad[i,:]).reshape(n_steps,1,nSHO)
    
            plt.plot(E_tij[0:200, 0, 0])
            plt.show()



if __name__=="__main__":
    main()
