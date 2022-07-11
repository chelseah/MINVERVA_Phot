#!/usr/bin/env python
import matplotlib 
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from setting import *
def create_reference_phot(target, target_mag, catalog, indir, dmag=0.5, ap='ap2_0', outlcpath=''):
    namelist = ["col%d" % (i+1) for i in range(4)]
    for i in range(4):
        for j in range(7):
            namelist+=["ap%d_%d" % (i, j)]

    usecols = ["col2", "col3", "col4"]
    for i in range(4):
        usecols+=["ap%d_%d" % (i,0)]
        usecols+=["ap%d_%d" % (i,1)]
        if i ==0:
            usecols+=["ap%d_%d" % (i,2)]
            usecols+=["ap%d_%d" % (i,3)]
        usecols+=["ap%d_%d" % (i,4)]
        usecols+=["ap%d_%d" % (i,5)]
        usecols+=["ap%d_%d" % (i,6)]
    
    flux_ref = []
    ticlist = []
    for index,row in catalog.iterrows():
        #print(row["mag"]-target_mag)
        tic = int(row["tic"])
        if tic == target:
            lcfile = os.path.join(indir, "%d.rlc" % tic)
            ticlist.append(tic)
        elif (row["mag"]>(target_mag-0.2)) and (row["mag"]<(target_mag+dmag)):
        #elif (row["mag"]>(target_mag-dmag)) and (row["mag"]<(target_mag+dmag)):
            #print(tic, row["mag"])
            lcfile = os.path.join(indir, "%d.rlc" % tic)
            #print(lcfile)
            if not os.path.exists(lcfile):
                continue
            ticlist.append(tic)
        else:
            continue
        
        df = pd.read_csv(lcfile, delim_whitespace=True, header=None, names=namelist, index_col=None, usecols=usecols)
        df = df.drop_duplicates()
        df = df.sort_values(by=['col2'])
        if flux_ref ==[]:
            
            flux_ref = 10**(np.array(df["%s_0" % ap])/-2.5)
        else:
            try:
                flux_ref += 10**(np.array(df["%s_0" % ap])/-2.5)
            except ValueError:
                continue
    #print (ticlist)
    datadic = {}
    for tic in ticlist:
        lcfile = os.path.join(indir, "%d.rlc" % tic)
        if not os.path.exists(lcfile):
            continue
        df = pd.read_csv(lcfile, delim_whitespace=True, header=None, names=namelist, index_col=None, usecols=usecols)
        df = df.drop_duplicates()
        df = df.sort_values(by=['col2'])
        cadence = np.array(df['col2'])
        flux = 10**(np.array(df["%s_0" % ap])/-2.5)
        if not (len(flux))==(len(flux_ref)):
            #print len(flux)
            continue
        flux/=(flux_ref-flux)
        flux/=np.nanmedian(flux)
        datadic["cadence"] = cadence
        datadic["%d" % tic] = flux
        x = df['%s_4' % ap]
        y = df['%s_5' % ap]
        datadic["x"] = x
        datadic["y"] = y
        print(tic, target)
        if tic == target:
            rms = np.nanmedian(np.abs(flux-np.nanmedian(flux)))
        #plt.plot(cadence, flux_ref/np.nanmedian(flux_ref), '.')
        #plt.plot(cadence, flux/np.nanmedian(flux), '.')
        #plt.show()
    fluxdic = pd.DataFrame.from_dict(datadic)
    fluxdic.to_csv(os.path.join(outlcpath, "%d.csv" % target), index=False)
    return [tic, rms] 


def select_clearence(target_ra, target_dec, target_mag, catalog, depth, maxsep=2.5/60.):
    c1 = SkyCoord(target_ra*u.deg, target_dec*u.deg, frame='icrs')
    dmag = np.log10(depth)*(-2.5)+0.1
    
    ticlist = []
    for index,row in catalog.iterrows():
        if row["mag"]-target_mag < dmag:
            c2 = SkyCoord(row["ra"]*u.deg, row["dec"]*u.deg, frame='icrs')
            sep = c1.separation(c2)
            if sep.deg<maxsep:
                stardic = {}
                stardic["tic"] = row["tic"]
                stardic["mag"] = row["mag"]
                stardic["sep"] = sep.deg
                ticlist.append(stardic)
    #print(ticlist)
    return ticlist


if __name__ == '__main__':
    #inlist = np.loadtxt(infile).astype(int)
    catfile = "RA%.1fDEC%.1f.cat" % (target_ra, target_dec)
    catalog = pd.read_csv(catfile, index_col=False, names=["tic", "ra", "dec", "mag", "col1", "col2", "col3", "col4", "col5"])
    #get reference photometry for the target
    indir = os.path.join(outpath, 'LC')
    
    create_reference_phot(target_id, target_mag, catalog, indir, dmag=refmagrange, ap=targetap, outlcpath=outlcpath)

    exit()
    tic_to_clear = select_clearence(target_ra, target_dec, target_mag, catalog, depth)
     
    for otherstars in tic_to_clear:
        tic, rms = create_reference_phot(otherstars["tic"], otherstars["mag"], catalog, indir, outlcpath=outlcpath)
        #print(tic,rms,otherstars["mag"],otherstars["sep"])

