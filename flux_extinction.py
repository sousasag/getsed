#!/usr/bin/python

##imports:

import matplotlib.pyplot as plt
import numpy as np
import getsed
from extinction import ccm89, apply

def flux_redden(wave, flux, av=1.0, rv=3.1):
    wavetot = wave.astype(np.double)
    fluxtot = flux.astype(np.double)
    flux_red = apply(ccm89(wavetot, av, rv), fluxtot)
    return flux_red






### Main program:
def main():

    wave, flux = getsed.get_sed_interpolated_cube(5777, 0.0, 4.4)
    fluxred = flux_redden(wave,flux)

    plt.plot(wave,flux)
    plt.plot(wave, fluxred)
    plt.show()

if __name__ == "__main__":
    main()
