from  astropy.io import fits
hdulist = fits.open("seds/kt06000g+4.5z+0.0sed.fits.gz", memmap=False)
hdulist = fits.open("sedskt06000g+4.5z+0.0_sed.fits.gz", memmap=False)
hdulist = fits.open("seds/kt06000g+4.5z+0.0_sed.fits.gz", memmap=False)
fits.info("seds/kt06000g+4.5z+0.0_sed.fits.gz", memmap=False)
hdulist[1].data
hdulist[3].data
%pylab
hdulist[2].data
hdulist[3].data
hdulist[3].data["TOTAL_FLUX"]
hdulist[3].data["TOTAL_FLUX"][0]
hdulist[1].data
hdulist[1].data["WAVELENGTH"]
hdulist[1].data["WAVELENGTH"][0]
plot(hdulist[1].data["WAVELENGTH"], hdulist[3].data["TOTAL_FLUX"])
plot(hdulist[1].data["WAVELENGTH"][0], hdulist[3].data["TOTAL_FLUX"])
wave = hdulist[1].data["WAVELENGTH"][0
]
wave.shape
wave = hdulist[1].data["WAVELENGTH"]
wave
wave.shape
flux = hdulist[3].data["TOTAL_FLUX"]
flux
flux.shape
flux = hdulist[3].data["TOTAL_FLUX"][0]
flux.shape
plot(wave, flux)
wave
plot(wave, flux)
xlim(4000,5000)
hdulist[1].columns[0].unit
hdulist[3].columns[0].unit
wave = hdulist[1].data["WAVELENGTH"] * 10000.
plot(wave, flux)
xlim(4000,5000)
from astropy import units as u
u.angstrom
wave = hdulist[1].data["WAVELENGTH"]
wave = hdulist[1].data["WAVELENGTH"] * u.micron
wave
wave.to(u.angstrom)
wave_a = wave.to(u.angstrom)
wave = hdulist[1].data["WAVELENGTH"] * u.micron
wave_a = wave.to(u.angstrom)
plot(wave_a, flux)
flux = hdulist[3].data["TOTAL_FLUX"][0] * u.mJy
flux_erg = flux.to(u.erg / u.cm **2. / u.s)
flux_erg = flux.to(u.erg / u.cm **2 / u.s)
nu = hdulist[1].data["FREQUENCY"]
hdulist[1].columns[1].unit
nu = hdulist[1].data["FREQUENCY"] * u.Hz
no
nu
hdulist[3].columns[0].unit
flux
flux.to(u.Jy)
flux = flux * nu
flux
flux = flux.to(u.erg / u.cm ** 2 / u.s)
flux
plot(wave_a, flux)
xlim(4000,5000)
flux
flux = hdulist[3].data["TOTAL_FLUX"][0] * u.mJy
flux
flux.to(u.Jy)
flux = flux.to(u.erg / u.cm ** 2 / u.s)
flux = flux * nu
flux
flux = flux.to(u.erg / u.cm ** 2 / u.s)
flux
flux = flux * wave_a
flux
flux = hdulist[3].data["TOTAL_FLUX"][0] * u.mJy
flux.to(u.Jy)
flux = flux * nu
flux
flux = flux.to(u.erg / u.cm ** 2 / u.s)
flux
flux = flux / wave_a
flux
flux = hdulist[3].data["TOTAL_FLUX"][0] * u.mJy
flux
flux = flux.to(u.erg / u.cm ** 2 / u.s / u.angstrom)
flux = flux.to(u.erg / u.cm ** 2 / u.s / u.Hz)
flux
flux = hdulist[3].data["TOTAL_FLUX"][0] * u.mJy
flux = flux.to(u.erg / u.cm ** 2 / u.s / u.Hz)
flux
%hist

