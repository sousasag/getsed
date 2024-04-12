#!/usr/bin/python

##imports:

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.table import Table
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import getbandpass as gbp
import getsed

## My functions:


def custumized_Simbad():
    customSimbad = Simbad()
    customSimbad.add_votable_fields('flux(G)') 
    customSimbad.add_votable_fields('flux(V)') 
    customSimbad.add_votable_fields('flux_error(G)') 
    customSimbad.add_votable_fields('flux_error(V)') 
    customSimbad.add_votable_fields('plx')
    customSimbad.add_votable_fields('plx_error')
    return customSimbad

def get_gaia_dr2_id(results_ids):
    for name in results_ids[::-1]:
        #print(name[0])
        if "Gaia DR2 " in name[0]:
            return name[0].split(" ")[-1]
    return -1

def get_gaiadr2(name):
    customSimbad=custumized_Simbad()
    if name[-2:] == " A":
        name =  name[:-2]
    if "(AB)" in name:
        name = name.replace("(AB)", "")
    if "Qatar" in name:
        name = name.replace("-","")
    try:
        result_ids = customSimbad.query_objectids(name)
    except:
        result_ids = customSimbad.query_objectids(name)
    if result_ids is None:
        gaiadr2 = -1
    else:
        gaiadr2 = get_gaia_dr2_id(result_ids)
    return gaiadr2

def get_2mass_data(gaiadr2, name_search=None, id2mass=None):
    vq2mass = Vizier(columns=['2MASS', 'RAJ2000', 'DEJ2000','Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag'], row_limit=5000)
    radius_search = 10.0*u.arcsec
    if name_search is None:
        result_2mass=vq2mass.query_object("Gaia DR2 "+str(gaiadr2), catalog=["II/246/out"], radius=radius_search)
    else:
        result_2mass=vq2mass.query_object(name_search, catalog=["II/246/out"], radius=radius_search)
    if id2mass is not None:
        il = np.where(result_2mass[0]['_2MASS'] == id2mass)[0][0]
        return result_2mass[0][il]
    return result_2mass

def get_allwise_data(gaiadr2, name_search=None, idallwise=None):
    vq2wise = Vizier(columns=['AllWISE', 'RAJ2000', 'DEJ2000', 'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag'], row_limit=5000)
    radius_search = 10.0*u.arcsec
    if name_search is None:
        result_allwise=vq2wise.query_object("Gaia DR2 "+str(gaiadr2), catalog=["II/328/allwise"], radius=radius_search)
    else:
        result_allwise=vq2wise.query_object(name_search, catalog=["II/328/allwise"], radius=radius_search)
    if idallwise is not None:
        il = np.where(result_allwise[0]['AllWISE'] == idallwise)[0][0]
        return result_allwise[0][il]
    return result_allwise

def get_gaiadr3_data(gaiadr2, gaia3, name_search=None):
    vq2 = Vizier(columns=['Source', '2MASS', 'AllWISE', 'RA_ICRS', 'DE_ICRS', 'Plx','e_Plx', 'FG','e_FG','Gmag','e_Gmag', 'BPmag','e_BPmag', 'RPmag','e_RPmag', 'o_Gmag'], row_limit=5000) 
    radius_search = 20.0*u.arcsec
    if name_search is None:
        result_gaia_vizier_dr3=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/355/gaiadr3"], radius=radius_search)
    else:
        result_gaia_vizier_dr3=vq2.query_object(name_search, catalog=["I/355/gaiadr3"], radius=radius_search*15.)
    try:
        iline3 = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaia3))[0][0]
        #print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline3])
    except IndexError:
        result_gaia_vizier_dr3=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/355/gaiadr3"], radius=radius_search*25.)
        try:
            iline3 = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaia3))[0][0]
        #    print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline3])
        except IndexError:
            result_gaia_vizier_dr3 =  Table( [[-1], [-1],[-1], [-1] ,[-1], [-1] ,[-1] ,[-1]   ,[-1] , [-1]    ] ,
                                    names=('Plx','e_Plx','Gmag','e_Gmag','FG','e_FG','RPmag','e_RPmag','BPmag','e_BPmag'))
            return result_gaia_vizier_dr3, -1
    #print("Iline2:", iline3)
    #print(result_gaia_vizier_dr3)
    return result_gaia_vizier_dr3[0][iline3]


def get_gaiaedr3_data(gaiadr2, gaia3, name_search=None):
    vq2 = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 'Plx','e_Plx', 'FG','e_FG','Gmag','e_Gmag', 'BPmag','e_BPmag', 'RPmag','e_RPmag', 'o_Gmag'], row_limit=5000) 
    radius_search = 20.0*u.arcsec
    if name_search is None:
        result_gaia_vizier_dr3=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/350/gaiaedr3"], radius=radius_search)
    else:
        result_gaia_vizier_dr3=vq2.query_object(name_search, catalog=["I/350/gaiaedr3"], radius=radius_search*15.)
    try:
        iline3 = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaia3))[0][0]
        #print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline3])
    except IndexError:
        result_gaia_vizier_dr3=vq2.query_object("Gaia DR2 "+str(gaiadr2), catalog=["I/350/gaiaedr3"], radius=radius_search*25.)
        try:
            iline3 = np.where(result_gaia_vizier_dr3[0]['Source'] == int(gaia3))[0][0]
        #    print(result_gaia_vizier_dr3[0]['Source','Plx','e_Plx', 'Gmag', 'e_Gmag', 'FG', 'e_FG', 'RPmag', 'e_RPmag', 'BPmag', 'e_BPmag'][iline3])
        except IndexError:
            result_gaia_vizier_dr3 =  Table( [[-1], [-1],[-1], [-1] ,[-1], [-1] ,[-1] ,[-1]   ,[-1] , [-1]    ] ,
                                    names=('Plx','e_Plx','Gmag','e_Gmag','FG','e_FG','RPmag','e_RPmag','BPmag','e_BPmag'))
            return result_gaia_vizier_dr3, -1
    #print("Iline2:", iline3)
    #print(result_gaia_vizier_dr3)
    return result_gaia_vizier_dr3[0][iline3]



def get_fluxEarth(gaiaid):

    data_gaia = get_gaiadr3_data(gaiaid, gaiaid)
    data_2mass = get_2mass_data(gaiaid, id2mass=data_gaia['_2MASS'])
    data_wise = get_allwise_data(gaiaid, idallwise=data_gaia['AllWISE'])


    bands = ["BP", "G", "RP", "J", "H", "K", "W1", "W2", "W3", "W4"]
    fluxes = []
    for b in bands:
        if b in ["BP", "G", "RP"]:
            data = data_gaia
        if b in ["J", "H", "K"]:
            data = data_2mass
        if b in ["W1", "W2", "W3", "W4"]:
            data = data_wise
        wave, fluxJy, fluxW, fluxerg = gbp.mag2flux(data[b+"mag"], b)
        print(b, wave, fluxJy, fluxW, fluxerg)
        fluxes.append((b, wave, fluxJy, fluxW, fluxerg))
    return fluxes


def get_fluxStar(teff, logg, feh):
    wave, flux = getsed.get_sed_interpolated_cube(teff, feh, logg)
    waveu, fluxu = getsed.get_sed_units(wave,flux)
    c = 29979245800.0 * u.cm / u.s
    flux = waveu**2/c * fluxu
    fluxW = flux.to(u.watt/u.m**2/u.Hz)
    fluxJy = flux.to(u.Jy)
    return waveu, fluxJy, fluxW, fluxu

def print_fluxes(gaiaid, teff, logg, feh):
    fluxes = get_fluxEarth(gaiaid)
    wave_s, fluxJy_s, fluxW_s, fluxu_s = get_fluxStar(teff, logg, feh)

    for b, w, fJym, fW, ferg in fluxes:
        fWs=fluxW_s[np.where(wave_s >= w.to(u.Angstrom))][0]
        theta = 2*np.sqrt(fW/fWs)*360*60*60/(2*np.pi)
        print (b, w, fJym, fW, ferg, fWs, theta)



### Main program:
def main():
    gaiadr2 = get_gaiadr2("HD80606")
    print(gaiadr2)
    print_fluxes(gaiadr2, 5543, 4.38, 0.31)


    return

#    data_gaia = get_gaiaedr3_data(gaiadr2, gaiadr2)
#    print(data_gaia)

    data_gaia = get_gaiadr3_data(gaiadr2, gaiadr2)
    print(data_gaia)

    print(data_gaia['_2MASS'])

    data_2mass = get_2mass_data(gaiadr2, id2mass=data_gaia['_2MASS'])
    print(data_2mass)

    print(data_gaia['AllWISE'])
    data_wise = get_allwise_data(gaiadr2, idallwise=data_gaia['AllWISE'])
    print(data_wise)

    fluxes_earth = get_fluxEarth(gaiadr2)




if __name__ == "__main__":
    main()
