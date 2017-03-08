#!/usr/bin/python
## Addapted from https://sedfitter.readthedocs.io

##imports:

from  astropy.io import fits
from astropy import units
import matplotlib.pyplot as plt



c = 29979245800.0 * units.cm / units.s
DISTANCE = 1. * units.kpc
PATH_SED = "/home/sousasag/Programas/GIT_projects/getsed/models_kurucz/seds/"

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


def read_ck04models(filename):
    """
    Physical fluxes of the spectra are given in FLAM surface flux units, 
    i.e. ergs cm^{-2} s^{-1} A^{-1}. These flux units differ from those in 
    the Castelli & Kurucz tables by a factor of 3.336 x 10^{-19} x lambda^{2} 
    x (4pi)^{-1}, i.e. are converted from ergs cm^{-2} s^{-1} Hz^{-1}steradian^{-1} 
    to ergs cm^{-2} s^{-1} A^{-1} by mutiplying the Castelli & Kurucz values by 
    3.336 x 10^{-19} x lambda^{2} x (4pi)^{-1}, where lambda is in Angstroms. To 
    convert to observed flux at Earth, multiply by a factor of (R/D)^2 where R is 
    the stellar radius, and D is the distance to Earth.
    """
    grav = filename.split('[')[-1].replace("]","")
    filein = filename.split('[')[0]
    data = fits.getdata(filein)
    print "reading: ", filein
    print "gravity: ", grav
    wave = data['WAVELENGTH'] * units.angstrom
    flux = data[grav] * units.erg / units.cm**2 / units.s / units.angstrom

    return wave, flux


def read_kurucz_sed(filename):
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

    flux_d = flux_c 

    return wave_a, flux_c





def plot_sed(wave, flux):
    plt.plot(wave, flux)
    plt.xlim(3000,10000)
    plt.show()


def oplotseds():
  test_sed = "models_kurucz/seds/kt06000g+4.5z+0.0_sed.fits.gz"
  wave,flux = read_kurucz_sed(test_sed)
  plt.plot(wave, flux)
  test_sed = "models_kurucz/seds/kt06500g+4.5z+0.0_sed.fits.gz"
  wave,flux = read_kurucz_sed(test_sed)
  plt.plot(wave, flux)
  test_sed = "models_kurucz/seds/kt07000g+4.5z+0.0_sed.fits.gz"
  wave,flux = read_kurucz_sed(test_sed)
  plt.plot(wave, flux)
  test_sed = "models_kurucz/seds/kt07500g+4.5z+0.0_sed.fits.gz"
  wave,flux = read_kurucz_sed(test_sed)
  plt.plot(wave, flux)
  test_sed = "models_kurucz/seds/kt08000g+4.5z+0.0_sed.fits.gz"
  wave,flux = read_kurucz_sed(test_sed)
  plt.plot(wave, flux)
  plt.xlim(3000,10000)
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

#  oplotseds()
#  return


  flux_Vega = 3836
  ll = 550 * units.nm
  change_flux(flux_Vega, ll)

  flux_test = 0.391483
  change_flux(flux_test, ll)

  test_sed = "models_kurucz/seds/kt07250g+4.5z+0.0_sed.fits.gz"
  test_sed = "models_kurucz/seds/kt08000g+2.5z-2.5_sed.fits.gz"
  test_sed = "models_kurucz/seds/kt04500g+4.0z+0.0_sed.fits.gz"

  wave,flux = read_kurucz_sed(test_sed)

  wave_n = wave.to(units.nm)
  flux_n = flux.to(units.erg / units.s / units.cm**2 / units.nm)

  plot_sed(wave, flux)


  test_sed = "ck04models/ckp00/ckp00_4500.fits[g40]"

  wave_hst , flux_hst = read_ck04models(test_sed)

  

  plot_sed(wave_hst, flux_hst)
  

if __name__ == "__main__":
    main()

