import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt

def main():
    nSHO = 357
    rand.seed(90211)
    #amp = np.array([np.exp(np.log(8) * abs(rand.normal())) for x in [0]*357])
    omega_cm = np.linspace(0,1500, 10000)
    beta = 1
    freq = 1 * np.tanh(beta * omega_cm / 2) *  ( ( 100 / (20000 + np.square(omega_cm)))  + \
                                                 ( 100 / (30000 + np.square(omega_cm)))  + \
                                                 ( 1000 / (2 * (10000 + np.square(omega_cm - 1000)))) + \
                                                 ( 1000 / (2 * (10000 + np.square(omega_cm - 500 ))))  )
    freq /= sum(freq)
    CDF = np.cumsum(freq)

    def inverse(omega,CDF,y):
        return omega[ np.argwhere(CDF >= y)[0] ]
        
    wSHO = [inverse(omega_cm, CDF, rand.rand()) for y in [0]*nSHO]
    H, bins = np.histogram(wSHO, bins=30)
    binsctr = .5 * bins[1:] + .5 * bins[:-1]

    plt.plot(binsctr, H/float(np.sum(H * (bins[1:] -  bins[:-1]))))
    plt.plot(omega_cm, freq/float(omega_cm[1] - omega_cm[0]))
    plt.show()



if __name__=="__main__":
    main()
