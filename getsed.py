#!/usr/bin/python
## some functions were addapted from https://sedfitter.readthedocs.io
## other grid taken from 
#    ftp://ftp.stsci.edu/cdbs/grid/ck04models/AA_README
#    ftp://ftp.stsci.edu/cdbs/grid/ck04models
#    http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html
#    http://www.stsci.edu/instruments/observatory/PDF/scs8.rev.pdf


DIR_SED = "/home/sousasag/Programas/GIT_projects/getsed/"

##imports:

from  astropy.io import fits
from astropy import units
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
from scipy.interpolate import RegularGridInterpolator

c = 29979245800.0 * units.cm / units.s
DISTANCE = 1. * units.kpc

NSEDS = 6688
NWAVE = 1221

## My functions:

def change_flux(flux_Jy, wave):
    """
    test with change_flux(3836, 550)

    https://www.gemini.edu/cgi-bin/sciops/instruments/michelle/magnitudes.pl?magnitude=9.978&wavelength=550&filter=Johnson+V&option=magnitude
    test with change_flux(3836, 550)
    """

    if type(flux_Jy) == units.quantity.Quantity:
        if not flux_Jy.unit == units.Unit("Jy"):
            flux_Jy = flux_Jy.to(units.Jy)
    else:
        flux_Jy = flux_Jy * units.Jy

    if type(wave) == units.quantity.Quantity:
        if not wave.unit == units.Unit("cm"):
            wave = wave.to(units.cm)
    else:
        wave = (wave * units.nm).to(units.cm)

    c = 29979245800.0 * units.cm / units.s

    print flux_Jy
    print wave
    print c

    Fn = flux_Jy.to(units.erg / units.s / units.cm**2 / units.Hz)
    Fl = Fn * c / wave**2.
    Fl = Fl.to(units.watt / units.m**2 / units.micron)
    print Fl
    Fl = Fl.to(units.erg / units.s / units.cm**2 / units.angstrom)
    lFl = Fl * wave
    lFl = lFl.to(units.watt / units.m **2)
    print Fl
    print lFl





def parse_units(unit_str):
    if unit_str == "MICRONS":
        return units.micron
    if unit_str == "mJy":
        return units.mJy
    if unit_str == "Jy":
        return units.Jy
    if unit_str == "HZ":
        return units.Hz
    if unit_str == "AU":
        return units.AU

def read_ck04models_numbers(filename):
    """
    Read a model from a fits file sed kurucz model

    """
    grav = filename.split('[')[-1].replace("]","")
    filein = filename.split('[')[0]
    data = fits.getdata(filein)
    print "reading: ", filein
    print "gravity: ", grav
    wave = data['WAVELENGTH']
    flux = data[grav]
    return wave, flux

def get_sed_units(wave,flux):
    """
    From the numbers read in ck04models and it adds units to the data
    """
    wave = wave * units.angstrom
    flux = flux * units.erg / units.cm**2 / units.s / units.angstrom
    return wave, flux    


def read_ck04models(filename):
    """
    Read a sed model in ck04models with units in the data (check Note)

    Note:

    Physical fluxes of the spectra are given in FLAM surface flux units, 
    i.e. ergs cm^{-2} s^{-1} A^{-1}. These flux units differ from those in 
    the Castelli & Kurucz tables by a factor of 3.336 x 10^{-19} x lambda^{2} 
    x (4pi)^{-1}, i.e. are converted from ergs cm^{-2} s^{-1} Hz^{-1}steradian^{-1} 
    to ergs cm^{-2} s^{-1} A^{-1} by mutiplying the Castelli & Kurucz values by 
    3.336 x 10^{-19} x lambda^{2} x (4pi)^{-1}, where lambda is in Angstroms. To 
    convert to observed flux at Earth, multiply by a factor of (R/D)^2 where R is 
    the stellar radius, and D is the distance to Earth.
    """
    wave, flux = read_ck04models_numbers(filename)
    return get_sed_units(wave, flux)


def read_kurucz_sed(filename):
    """
    Read a sed model in kurucz models (as in sed fitter)

    """

    hdulist = fits.open(filename, memmap=False)

    wave_unit = hdulist[1].columns[0].unit
    nu_unit   = hdulist[1].columns[1].unit
    flux_unit = hdulist[3].columns[0].unit
    ap_unit   = hdulist[2].columns[0].unit
    flux_er_unit = hdulist[3].columns[1].unit

    wave = hdulist[1].data.field("WAVELENGTH") * parse_units(wave_unit)
    nu = hdulist[1].data.field("FREQUENCY") * parse_units(nu_unit)
    ap = hdulist[2].data.field("APERTURE") * parse_units(ap_unit)
    flux = hdulist[3].data.field("TOTAL_FLUX")[0] * parse_units(flux_unit)
    flux_er = hdulist[3].data.field("TOTAL_FLUX_ERR")[0] * parse_units(flux_er_unit)

    #conversion of units to Angstrom and erg sec-1 cm-2 A-1
    flux = flux.to(units.Jy)
    flux_c = flux * c / wave.to(units.cm)**2
    flux_c = flux_c.to(units.erg / units.s / units.cm**2 / units.angstrom)

    wave_a = wave.to(units.angstrom)
    return wave_a, flux_c


def plot_sed(wave, flux):
    """
    simple plot of the sed model in range 3000-10000 Angstrom
    """
    plt.plot(wave, flux)
    plt.xlim(3000,10000)
    plt.show()


def oplotseds():
    """
    Replicate (out scalled) Fig 2.  Effective Temperature Determination (Niemczura book wroclaw)
    """
    test_sed = DIR_SED + "models_kurucz/seds/kt06000g+4.5z+0.0_sed.fits.gz"
    wave,flux = read_kurucz_sed(test_sed)
    plt.plot(wave, flux)
    test_sed = DIR_SED + "models_kurucz/seds/kt06500g+4.5z+0.0_sed.fits.gz"
    wave,flux = read_kurucz_sed(test_sed)
    plt.plot(wave, flux)
    test_sed = DIR_SED + "models_kurucz/seds/kt07000g+4.5z+0.0_sed.fits.gz"
    wave,flux = read_kurucz_sed(test_sed)
    plt.plot(wave, flux)
    test_sed = DIR_SED + "models_kurucz/seds/kt07500g+4.5z+0.0_sed.fits.gz"
    wave,flux = read_kurucz_sed(test_sed)
    plt.plot(wave, flux)
    test_sed = DIR_SED + "models_kurucz/seds/kt08000g+4.5z+0.0_sed.fits.gz"
    wave,flux = read_kurucz_sed(test_sed)
    plt.plot(wave, flux)
    plt.xlim(3000,10000)
    plt.show()


def oplotseds2():
    """
    Replicate Fig 2.  Effective Temperature Determination (Niemczura book wroclaw)
    """
    test_sed = DIR_SED + "ck04models/ckp00/ckp00_6000.fits[g45]"
    wave,flux = read_ck04models(test_sed)
    plt.plot(wave, flux)
    test_sed = DIR_SED + "ck04models/ckp00/ckp00_6500.fits[g45]"
    wave,flux = read_ck04models(test_sed)
    plt.plot(wave, flux)
    test_sed = DIR_SED + "ck04models/ckp00/ckp00_7000.fits[g45]"
    wave,flux = read_ck04models(test_sed)
    plt.plot(wave, flux)
    test_sed = DIR_SED + "ck04models/ckp00/ckp00_7500.fits[g45]"
    wave,flux = read_ck04models(test_sed)
    plt.plot(wave, flux)
    test_sed = DIR_SED + "ck04models/ckp00/ckp00_8000.fits[g45]"
    wave,flux = read_ck04models(test_sed)
    plt.plot(wave, flux)
    plt.xlim(3000,10000)
    plt.show()


def read_save_grid(filebin = 'sed_data.pkl', val_type = 'float32'):
    """
    Read from fits and writting the fill sed grid into a single binary file
    """
    data_grid = fits.getdata(DIR_SED + "ck04models/catalog.fits")
    sed_grid = np.zeros((NSEDS*NWAVE,5))
    t0 = time.time()
    i=0
    count = 0
    for line in data_grid:
        teff, met, logg = line['INDEX'].split(',')
        teff = float(teff)
        met  = float(met)
        logg = float(logg)
        wave, flux = read_ck04models_numbers(DIR_SED + 'ck04models/'+line['FILENAME'])
        print i, len(data_grid),line['INDEX'], line['FILENAME'], teff, met, logg,len(wave), wave[0]
        if np.all(flux == 0):
            print "skipping zero flux sed:"
        else:
            for j in range(len(wave)):
                #print line['FILENAME'], teff, met, logg, wave[j], flux[j]
                sed_grid[count,:] = [teff,met,logg,wave[j],flux[j]]
                count +=1
        i += 1
    sed_grid = sed_grid[:count,:].astype(val_type)

    t1 = time.time()
    total2 = t1-t0
    print total2

    t0 = time.time()
    output = open('sed_data.pkl', 'wb')
    pickle.dump(sed_grid, output,-1)
    output.close()
    t1 = time.time()
    total2 = t1-t0
    print total2
    return sed_grid


def read_grid_pickle(filebin = 'sed_data.pkl'):
    """
    read a binary file with the grid created by the read_save_grid function
    """
    t0 = time.time()
    inputfile = open(filebin, 'rb')
    sed_grid = pickle.load(inputfile)
    inputfile.close()
    t1 = time.time()
    total3 = t1-t0
    print total3
    return sed_grid




def get_sed_interpolated_cube(teff, met, logg):
    """
    Returns an interpolated sed model:

    Args:
        teff: effective temperature
        met: metallicity
        logg: surface gravity

    Return:
        wave, flux - tuple with numpy array with the interpolated sed

    Examples:
        >>> wave, flux = get_sed_interpolated_cube(5777, 0, 4.44)
    """

    data_grid = fits.getdata(DIR_SED + "ck04models/catalog.fits")
    teffv = []
    metv  = []
    loggv = []
    for line in data_grid:
        teffi, meti, loggi = line['INDEX'].split(',')
        teffv.append(float(teffi))
        metv.append(float(meti))
        loggv.append(float(loggi))
    teff_u = np.unique(teffv)
    met_u = np.unique(metv)
    logg_u = np.unique(loggv)

# Getting the points of the cube to interpolate
    teff_l = teff_u[np.where(teff_u<teff)[0][-1]]
    teff_h = teff_u[np.where(teff_u>=teff)[0][0]]
    met_l = met_u[np.where(met_u<met)[0][-1]]
    met_h = met_u[np.where(met_u>=met)[0][0]]
    logg_l = logg_u[np.where(logg_u<logg)[0][-1]]
    logg_h = logg_u[np.where(logg_u>=logg)[0][0]]
    print teff_l, teff, teff_h
    print met_l, met, met_h
    print logg_l, logg, logg_h

    list_cube = [(teff_l, met_l, logg_l),(teff_l, met_l, logg_h),
                 (teff_l, met_h, logg_l),(teff_l, met_h, logg_h),
                 (teff_h, met_l, logg_l),(teff_h, met_l, logg_h),
                 (teff_h, met_h, logg_l),(teff_h, met_h, logg_h)]

# Reading the seds in the cube
    list_sed = []
    for t,m,l in list_cube:
        if m >=0:
            sm = "p%02d" % int(abs(m)*10.)
        else:
            sm = "m%02d" % int(abs(m)*10.)
        sl = "%2d" % int(l*10.)
        if t > 9999:
            st = "%5d" % int(t)
        else:
            st = "%4d" % int(t)
        file_name = DIR_SED + "ck04models/ck"+sm+"/ck"+sm+"_"+st+".fits[g"+sl+"]"
        wave, flux = read_ck04models_numbers(file_name)
        if np.all(flux == 0):
            print "Problem with sed: ", file_name
            raise ValueError('Sed in interpolation cube with zero flux values')
        list_sed.append((wave,flux))


# Interpolating the sed
    wave_i = []
    flux_i = []
    t = np.linspace(teff_l, teff_h, 2)
    m = np.linspace(met_l, met_h, 2)
    l = np.linspace(logg_l, logg_h, 2)
    V = np.zeros((2,2,2))
    pt = (teff, met, logg)

    for i in range(NWAVE):
        V[0,0,0] = list_sed[0][1][i]
        V[0,0,1] = list_sed[1][1][i]
        V[0,1,0] = list_sed[2][1][i]
        V[0,1,1] = list_sed[3][1][i]
        V[1,0,0] = list_sed[4][1][i]
        V[1,0,1] = list_sed[5][1][i]
        V[1,1,0] = list_sed[6][1][i]
        V[1,1,1] = list_sed[7][1][i]
        fn = RegularGridInterpolator((t,m,l), V)
        flux_i.append(fn(pt))
        wave_i.append(list_sed[0][0][i])
    wave_i = np.array(wave_i)
    flux_i = np.array(flux_i)

    return wave_i, flux_i



def test_interpolation():
    """
    Compare this interpolation with a 
    previous interpolation generated with iuerdaf kurget
    """
    wave_i, flux_i = get_sed_interpolated_cube(8075, 0.35, 4.9)
    wave_t, flux_t = np.loadtxt(DIR_SED + 'kuruczbm.dat', unpack = True)
    scale = np.mean(flux_i)/np.mean(flux_t)
#    plt.plot(wave_i, flux_i, linewidth=3, color='k')
#    plt.plot(wave_t, flux_t * scale, linewidth=3, color='g')
    plt.plot(wave_i, (flux_i - flux_t*scale)/flux_i)
    plt.xlim(3000,11000)
    plt.show()


def compare_grids():
    """
    compare the 2 format grids
    """
    test_sed = DIR_SED + "models_kurucz/seds/kt08000g+2.5z-2.5_sed.fits.gz"

    wave_1,flux_1 = read_kurucz_sed(test_sed)
    wave_1 = wave_1[:NWAVE]
    flux_1 = flux_1[:NWAVE]

    test_sed = DIR_SED + "ck04models/ckm25/ckm25_8000.fits[g25]"

    wave_2, flux_2 = read_ck04models(test_sed)

    scale = np.mean(flux_2)/np.mean(flux_1)



#    print flux_1.shape, flux_2.shape
#    print (flux_1 - flux_2)[300:600]
#    plt.plot(wave_2, flux_2, linewidth=3, color='k')
#    plt.plot(wave_1, flux_1 * scale, linewidth=3, color='g')
    plt.plot(wave_1, (flux_2 - flux_1 * scale)/flux_2)
    plt.xlim(3000,11000)
    plt.show()


### Main program:
def main():
    """
    ftp://ftp.stsci.edu/cdbs/grid/ck04models/AA_README
    ftp://ftp.stsci.edu/cdbs/grid/ck04models
    http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html
    http://www.stsci.edu/instruments/observatory/PDF/scs8.rev.pdf
    """
    print "Hello"

#    compare_grids()
#    return
#    test_interpolation()
#    return
#    oplotseds()
#    return
#    oplotseds2()
#    return

  

if __name__ == "__main__":
    main()

