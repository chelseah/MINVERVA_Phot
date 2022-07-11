#!/usr/bin/env python
from __future__ import print_function
import os
import sys
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '..'))
import astropy
import astropy.wcs as wcs
import phot
from phot import Sextractor,Fistar,Anet,Starlist
import pandas as pd
from setting import *
import sys
import numpy as np

def main(infile,outdir,ra0,dec0,astrom_dir,cfgfile,cam):
    print("doing photometry on file %s" % infile)
    #infile = os.path.join(inpath, infile)
    #source_extractor =  Sextractor()
    


    if cam == "QHY600":
        threshold = 30
        width = 9600
        height = 6422
        radius= 0.3
        upper = 0.2
        lower = 0.1

    elif cam == "ZWO295":
        threshold =10
        width = 4144
        height = 2822
        radius = 1.
        upper = 0.3
        lower = 0.1

    elif cam == "andor":
        print("camera is andor!!!!!!!!!!!!!!!!!!!!!!")
        threshold = 100
        width = 2048
        height = 2048
        radius = 0.545
        upper = 1.0
        lower = 0.6
    else:
        raise NotImplementedError("Camera type not found! Call Chelsea")

    #@width=5424
    #height= 3638
    #radius = 0.3
    source_extractor =  Sextractor(threshold=threshold)
    source_extractor(infile, out_dir=astrom_dir)
    starlist = Starlist(os.path.join(astrom_dir, source_extractor.get_outfile(infile)), colx=4, coly=5, colmag=2) 
    astrom = Anet(ra=ra0, dec=dec0, q=0.1, tweak=2, order=4, upper=upper,lower=lower, radius=radius, catfile=None, maxdistance=1, width=width, height=height)
    
    astrom(starlist, out_dir=outdir, continuesolve=False, refine=False, source="sex", astromcfg=cfgfile)
    wcsfile = os.path.join(astrom_dir, astrom.get_wcsfile(os.path.join(astrom_dir,starlist.name)))
    #break
    print(astrom_dir, wcsfile)

    w = wcs.WCS(wcsfile)
    wcsheader = w.to_header()
    try:
        costheta = wcsheader['PC1_2']
        sintheta = wcsheader['PC1_1']
        theta = np.arctan2(costheta,sintheta)
        theta = 180 - theta*180/np.pi
    except KeyError:
        theta = np.nan
    ra0, dec0 = w.all_pix2world(width/2.,height/2., 0)



    
    print(w)
    print(ra0,dec0)

    return astrom_dir,wcsfile,w,ra0,dec0,theta


if __name__=='__main__':
    infile = sys.argv[1]
    ra0 = float(sys.argv[2])
    dec0 = float(sys.argv[3])
    astrom_dir = sys.argv[4]
    cfgfile = sys.argv[5]
    cam = sys.argv[6]
    main(infile,ra0,dec0,astrom_dir,cfgfile,cam)
    
