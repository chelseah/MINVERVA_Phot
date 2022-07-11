#!/usr/bin/env python
from barycorrpy import utc_tdb
import astropy
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
import pandas as pd
import numpy as np
import os
from setting import *
df = pd.read_csv(inlist, index_col=False)


star = SkyCoord(target_ra*u.deg, target_dec*u.deg, frame='icrs')
df["airmass"] = np.zeros(len(df.fitsfile))
df["jd"] = np.zeros(len(df.fitsfile))
df["exp"] = np.zeros(len(df.fitsfile))
df["hjd"] = np.zeros(len(df.fitsfile))
df["bjd"] = np.zeros(len(df.fitsfile))
df["temp"] = np.zeros(len(df.fitsfile))
df["filter"] = np.empty(len(df.fitsfile), dtype='S256')

for index,row in df.iterrows():
    infile = os.path.join(inpath, row["fitsfile"])
    hdulist = fits.open(infile)
    df.at[index, "jd"] = hdulist[0].header["JD"] 
    jd = float(hdulist[0].header["JD"])
    JDUTC = Time(jd, format='jd', scale='utc')
    try:
        df.at[index, "exp"] = hdulist[0].header["EXPOSURE"] 
    except KeyError:
        pass
    try:
        df.at[index, "hjd"] = hdulist[0].header["JD-HELIO"]
    except KeyError:
        pass
    try:
        df.at[index, "filter"] = hdulist[0].header["FILTER"]
    except KeyError:
        pass
    try:
        df.at[index, "temp"] = hdulist[0].header["CCD-TEMP"]
    except KeyError:
        pass
    sitelat = -27.7977
    sitelong = 151.8554
    siteheight = 682
    site = EarthLocation(lat=sitelat*u.deg, lon=sitelong*u.deg, height=siteheight*u.m)
    altz = star.transform_to(AltAz(obstime=JDUTC, location=site))
    airmass = altz.secz
    df.at[index,"airmass"] = airmass 
    bjdtdb = utc_tdb.JDUTC_to_BJDTDB(JDUTC,ra=target_ra ,dec=target_dec, lat=sitelat, longi=sitelong, alt=siteheight, ephemeris='de405.bsp')[0][0]
    df.at[index, "bjd"] = bjdtdb


df.to_csv("fitsheader.csv", index=False)
