#!/usr/bin/python

##imports:

import matplotlib.pyplot as plt
import numpy as np


## My functions:

file_bp_gaiadr3 = "bandpasses/GaiaEDR3_passbands_zeropoints_version2/passband.dat"

dype_gdr3 = [('lambda', '<f8'),('GPb', '<f8')  , ('e_GPb', '<f8')  , 
                               ('BPPb', '<f8') , ('e_BPPb', '<f8') , 
                               ('RPPb', '<f8') , ('e_RPPb', '<f8')  ]

def get_passband_gaia_dr3():
    data = np.loadtxt(file_bp_gaiadr3, dtype = dtype_gdr3)
    for col,_ in dtype_gdr3:
        data[col][data[col] == 99.99] = 0
    return data


### Main program:
def main():
    print ("Hello")



if __name__ == "__main__":
    main()
