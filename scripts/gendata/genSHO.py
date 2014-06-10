import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import argparse
import h5py

def SHOparam(s):
    try:
        tau, Asq = map(float,s.split(','))
        return tau,np.sqrt(Asq)
    except:
        raise argparse.ArgumentTypeError("Argument must be tau,Asq")

def main():
    parser = argparse.ArgumentParser(description = "Generate SHO spectral density bath data to be used as proxy data for understanding the results of pyPCA module analysis tools. [p(A) propto (amppeak / A^pow)]")
    parser.add_argument('-nSHO', type=int, default=357, help="Number of harmonic oscillators to use in the dynamics")
    parser.add_argument('-dt', type=float, default=.01, help="Timestep for simulation, ps. Stored in hdf5 file.")
    parser.add_argument('-T', type=float, default=2000., help="Total time to generate, ps.")
    parser.add_argument('-h5out', help="Location for the output of the dynamics")
    parser.add_argument('-h5dset', default='alltimes', help="Name of the timeseries in the hdf5 database")
    parser.add_argument('-DOF', type=int, default=1, help="Number of DOF merged into each time slice")
    parser.add_argument('-pow', type=float, default=1., help="Power law distribution for the amplitudes [p(A) propto (amppeak / A^pow)]")
    parser.add_argument('-amppeak', type=float, default=100., help="Peak amplitude (i.e. the scale of amplitude decay) for amplitude distribution [p(A) propto (amppeak / A^pow)]")
    parser.add_argument('-dynamics', choices=['microcanonical','brownian'], default='microcanonical', help="Set the type of dynamics to generate your data with")
    parser.add_argument('-setSHO', type=SHOparam, nargs="*", help="Set SHO parameters explicitly.")
    args = parser.parse_args()

    if args.dynamics == 'brownian' and args.DOF!=1:
        raise NotImplementedError("WARNING: DOF (={}) >1 and brownian dynamics set. Brownian dynamics are only enabled to use one dof. Please consider a patch.".format(args.DOF))

    n_steps = round(args.T / args.dt) + 1
    
    if len(args.setSHO) != 0:
        nSHO = len(args.setSHO)
        args.setSHO = np.array(args.setSHO)
        aSHO = args.setSHO[:,1]
        wSHO_ps = np.reciprocal(args.setSHO[:,0])
        wSHO_ps = wSHO_ps.reshape([1,nSHO])
        plt.plot(np.reciprocal(wSHO_ps[0,:]))
        plt.ylabel("tau, ps")
        plt.xlabel("mode number")
        plt.show()
        #raise NotImplementedError("Did not finish implementation of manually set SHO parameters.")
    else:
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
        A  = np.linspace(0.,args.amppeak, args.amppeak+1)
        pA = args.amppeak * np.reciprocal(np.power(A,args.pow))
        pA[0] = 0
        CDF_A  = np.cumsum(pA / np.sum(pA))
        #aSHO     = rand.gamma(10., .1, (args.DOF, nSHO))
        aSHO     = np.interp(rand.rand(nSHO), CDF_A, A)
        print max(aSHO.flatten())
        H, bins = np.histogram(aSHO, bins=400)
        binsctr = .5 * bins[1:] + .5 * bins[:-1]
        plt.plot(binsctr, H/float(np.sum(H * (bins[1:] -  bins[:-1]))))
        if (args.pow > 1):
            C = args.amppeak
            norm = .5 * C + (args.pow-1) * (C - C**(2. - args.pow))
            plt.plot(A, pA/norm)
        plt.ylabel("p(A)")
        plt.xlabel("Amplitude")
        plt.show()
        wSHO_ps  = np.interp(rand.rand(args.DOF, nSHO), CDF_W, omega_cm/planck_cmps) 
    
    pSHO_rad = rand.rand(wSHO_ps.shape[0],wSHO_ps.shape[1]) *2.*np.pi
    print "Shapes for A, Phase, and Freq:", aSHO.shape, pSHO_rad.shape, wSHO_ps.shape

    H, bins = np.histogram(wSHO_ps, weights=aSHO + np.zeros(wSHO_ps.shape), bins=50)
    binsctr = .5 * bins[1:] + .5 * bins[:-1]
    plt.plot(binsctr, H/float(np.sum(H * (bins[1:] -  bins[:-1]))))
    if len(args.setSHO) == 0:
        plt.plot(omega_cm/planck_cmps, freq / float(omega_cm[1] - omega_cm[0])*planck_cmps)
    plt.xlabel("1/ps, frequency")
    plt.ylabel("Density of SHO")
    plt.title("Spectral density vs actual modes")

    plt.show()

    t = np.linspace(0,args.T, n_steps)

    if args.h5out and args.dynamics=='microcanonical':
        with h5py.File(args.h5out,'w') as f:
            E_tij = f.create_dataset(args.h5dset, (n_steps, 1, nSHO), fillvalue=0.)
            E_tij.attrs['dt'] = args.dt
            for i in xrange(args.DOF):
                E_tij[:] += (aSHO * np.cos(np.outer(t,wSHO_ps[i,:]) + pSHO_rad[i,:])).reshape(n_steps,1,nSHO)
    
            plt.plot(E_tij[0:200, 0, 0])
            plt.show()

    if args.h5out and args.dynamics=='brownian':
        with h5py.File(args.h5out,'w') as f:
            print("Running dynamics!")
            E_tij = f.create_dataset(args.h5dset, (n_steps, 1, nSHO), fillvalue=0.)
            E_tij.attrs['dt'] = args.dt
            x_j  = rand.normal(size=(1,nSHO))
            plt.plot(x_j[0,:])
            plt.show()
            dt_j = args.dt * wSHO_ps
            print wSHO_ps
            E_tij[0,0,:] = x_j * aSHO
            print "xj, dtj, nSHO",x_j.shape, dt_j.shape, nSHO
            for t in np.arange(1,n_steps, dtype=int):
                if t%500==0:
                    print "t={} of {}".format(t,n_steps)
                x_j = x_j * np.exp(-dt_j) + rand.normal(size=(1,nSHO)) * np.sqrt(1 - np.exp(-2 * dt_j))
                E_tij[t,0,:] = x_j * aSHO

            plt.plot(np.linspace(0,args.T,n_steps)[0:2000],E_tij[0:2000, 0, 0])
            plt.show()



if __name__=="__main__":
    main()
