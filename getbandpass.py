#!/usr/bin/python

##imports:

import matplotlib.pyplot as plt
import numpy as np

#JHK
#https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
## My functions:

file_bp_gaiadr3 = "bandpasses/GaiaEDR3_passbands_zeropoints_version2/passband.dat"
file_2mass_J = "bandpasses/2MASS_J.dat"
file_2mass_H = "bandpasses/2MASS_H.dat"
file_2mass_Ks = "bandpasses/2MASS_Ks.dat"

file_wises = ['bandpasses/RSR-W1.txt', 'bandpasses/RSR-W2.txt', 'bandpasses/RSR-W3.txt', 'bandpasses/RSR-W4.txt']

dtype_gdr3 = [('lambda', '<f8'),('GPb', '<f8')  , ('e_GPb', '<f8')  , 
                               ('BPPb', '<f8') , ('e_BPPb', '<f8') , 
                               ('RPPb', '<f8') , ('e_RPPb', '<f8')  ]



def get_passband_gaia_dr3():
    data = np.loadtxt(file_bp_gaiadr3, dtype = dtype_gdr3)
    for col,_ in dtype_gdr3:
        data[col][data[col] == 99.99] = 0
    return data

def get_passband_2MASS():
    J = np.loadtxt(file_2mass_J, skiprows=1, dtype = [('lambda', '<f8'), ('TP', '<f8')])
    H = np.loadtxt(file_2mass_H, skiprows=1, dtype = [('lambda', '<f8'), ('TP', '<f8')])
    Ks = np.loadtxt(file_2mass_Ks, skiprows=1, dtype = [('lambda', '<f8'), ('TP', '<f8')])
    return J, H, Ks

def get_passband_WISE():
    Wi_list = []
    for f in file_wises:
        Wi = np.loadtxt(f, skiprows=1, usecols=(0,1), dtype = [('lambda', '<f8'), ('TP', '<f8')])
        Wi_list.append(Wi)
    return tuple(Wi_list)


#Table 1 - 2MASS Isophotal Bandpasses and Fluxes-for-0-magnitude from Cohen et al. (2003)
#Band    Lambda (µm) Bandwidth (µm)  Fnu - 0 mag (Jy)    Flambda - 0 mag (W cm-2 µm-1)
#J   1.235 ± 0.006   0.162 ± 0.001   1594  ± 27.8    3.129E-13 ± 5.464E-15
#H   1.662 ± 0.009   0.251 ± 0.002   1024  ± 20.0    1.133E-13 ± 2.212E-15
#Ks  2.159 ± 0.011   0.262 ± 0.002   666.7 ± 12.6    4.283E-14 ± 8.053E-16

def mag2flux(mag, band):
    if band == "Ks":
        wave = 2159
        zero = 667
    flux = zero * 10**(mag/2.5)
    return flux

### Main program:
def main():

    data_gaia = get_passband_gaia_dr3()
    J, H, Ks = get_passband_2MASS()
    W1, W2, W3, W4 = get_passband_WISE()


    fig_gaia = plt.figure()
    ax1 = fig_gaia.add_subplot(311)
    ax2 = fig_gaia.add_subplot(312)
    ax3 = fig_gaia.add_subplot(313)
    ax1.plot(data_gaia['lambda'], data_gaia["GPb"])
    ax2.plot(data_gaia['lambda'], data_gaia["BPPb"])
    ax3.plot(data_gaia['lambda'], data_gaia["RPPb"])
    plt.show()



    fig_2MASS = plt.figure()
    ax21 = fig_2MASS.add_subplot(311)
    ax22 = fig_2MASS.add_subplot(312)
    ax23 = fig_2MASS.add_subplot(313)
    ax21.plot(J['lambda'], J['TP'])
    ax22.plot(H['lambda'], H['TP']) 
    ax23.plot(Ks['lambda'], Ks['TP'])
    plt.show()

    fig_WISE = plt.figure()
    ax31 = fig_WISE.add_subplot(111)
    ax31.plot(W1['lambda'], W1['TP'])
    ax31.plot(W2['lambda'], W2['TP'])
    ax31.plot(W3['lambda'], W3['TP'])
    ax31.plot(W4['lambda'], W4['TP'])
    ax31.set_xlim(2,30)
    ax31.set_ylim(0,1)
    plt.show()





if __name__ == "__main__":
    main()
