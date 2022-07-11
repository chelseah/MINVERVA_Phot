#!/usr/bin/env python
from __future__ import print_function
import os
import sys
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '..'))
import astropy
import astropy.wcs as wcs
import phot
from phot import Catalog,Fistar,Anet,Fiphot, Starlist
from setting_astrom import *

if __name__=='__main__':
    astrom_dir = os.path.join(outpath,"astrom/")
    infile = os.path.join(inpath, sys.argv[1])
    source_extractor = Fistar(threshold=1000, sdkflag=True, fiign=False)
    source_extractor(infile, out_dir=astrom_dir)
    starlist = Starlist(os.path.join(astrom_dir, source_extractor.get_outfile(infile)), colx=2, coly=3, colmag=9) 
    astrom = Anet(ra=target_ra, dec=target_dec, q=0.1, tweak=2, order=4, radius=radius, catfile=None, maxdistance=1, width=width, height=height)
    
    astrom(starlist, out_dir=astrom_dir, continuesolve=False, refine=False)
    wcsfile = os.path.join(astrom_dir, astrom.get_wcsfile(os.path.join(astrom_dir,starlist.name)))
    #break
    if not os.path.exists(wcsfile):
        exit()
    print(astrom_dir, wcsfile)
    w = wcs.WCS(wcsfile)
    ra0, dec0 = w.all_pix2world(width/2.,height/2., 0)
    print(ra0,dec0) 
