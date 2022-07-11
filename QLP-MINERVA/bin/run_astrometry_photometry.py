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

DEBUG=True

def main(infile,outdir,ra0,dec0,astrom_dir,cfgfile,fitshape):
    
    if DEBUG:
        print("doing photometry on file %s" % infile)
    #infile = os.path.join(inpath, infile)
    source_extractor =  Sextractor()
    source_extractor =  Sextractor(threshold=10)
    # TBD: need to change to what works with the raw frames
    #threshold = 1000
    #source_extractor = Fistar(threshold = threshold, sdkflag = True, fiign = False)
    source_extractor(infile, out_dir=astrom_dir)
    width=fitshape[1]
    height=fitshape[0]

    radius=0.2
    if max(fitshape) > 9000:
        radius = 0.27
    if max(fitshape) < 5000:
        #TBD:need to verify the number is correct
        radius = 0.115

    starlist = Starlist(os.path.join(astrom_dir, source_extractor.get_outfile(infile)), colx=4, coly=5, colmag=2) 
    #starlist = Starlist(os.path.join(astrom_dir, source_extractor.get_outfile(infile)), colx=2, coly=3, colmag=9) 
    astrom = Anet(ra=ra0, dec=dec0, q=0.1, tweak=2, order=4, radius=radius, catfile=None, maxdistance=1, width=width, height=height)
    
    astrom(starlist, out_dir=outdir, continuesolve=False, refine=False, source="sex", astromcfg=cfgfile)
    wcsfile = os.path.join(astrom_dir, astrom.get_wcsfile(os.path.join(astrom_dir,starlist.name)))
    #break
    #print("astrom_dir,wcs_file",astrom_dir, wcsfile)
    w = wcs.WCS(wcsfile)
    wcsheader = w.to_header()
    costheta = wcsheader['PC1_2']
    sintheta = wcsheader['PC1_1']
    theta = np.arctan2(costheta,sintheta)
    theta = 180 - theta*180/np.pi
    ra0, dec0 = w.all_pix2world(width/2.,height/2., 0)
    #print("w",w)
    #print("ra0,dec0",ra0,dec0)

    return astrom_dir,wcsfile,w,ra0,dec0,theta


if __name__=='__main__':
    infile = sys.argv[1]
    ra0 = float(sys.argv[2])
    dec0 = float(sys.argv[3])
    astrom_dir = sys.argv[4]
    cfgfile = sys.argv[5]
    main(infile,ra0,dec0,astrom_dir,cfgfile)
    
