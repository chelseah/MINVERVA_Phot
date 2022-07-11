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
import pandas as pd
from setting import *

if __name__=='__main__':
    df = pd.read_csv(inlist, index_col=False)
    astrom_dir = os.path.join(outpath,"astrom/")
    for index,row in df.iterrows():
        print("doing photometry on file %s" % row)
        infile = os.path.join(inpath, row["fitsfile"])
        framenb = row["frameno"] 
        source_extractor = Fistar(threshold=1000, sdkflag=True, fiign=False)
        source_extractor(infile, out_dir=astrom_dir)
        starlist = Starlist(os.path.join(astrom_dir, source_extractor.get_outfile(infile)), colx=2, coly=3, colmag=9) 
        astrom = Anet(ra=target_ra, dec=target_dec, q=0.1, tweak=2, order=4, radius=radius, catfile=None, maxdistance=1, width=width, height=height)
        
        astrom(starlist, out_dir=astrom_dir, continuesolve=False, refine=False)
        wcsfile = os.path.join(astrom_dir, astrom.get_wcsfile(os.path.join(astrom_dir,starlist.name)))
        #break
        if not os.path.exists(wcsfile):
            continue
        print(astrom_dir, wcsfile)
        w = wcs.WCS(wcsfile)
        ra0, dec0 = w.all_pix2world(width/2.,height/2., 0)
        print(ra0,dec0)
        catfile = "RA%.1fDEC%.1f.cat" % (ra0, dec0)
        xylist = os.path.join(astrom_dir,os.path.splitext(row["fitsfile"])[0]+".xyls")
        gaiastars = Catalog(catfile, ra0, dec0, 25/60., colid=1, colra=2, coldec=3, colmag=4, colx=5, coly=6, method="TIC")
        if not os.path.exists(catfile):
            gaiastars.query()
        #astrom.proj_sky_to_xy(xylist, inputcat=catfile,xlim=[0,2048], ylim=[0,2048])
        astrom.proj_sky_to_xy(xylist, inputcat=catfile,xlim=[0, width], ylim=[0,height])
        catproj = Starlist(xylist)
        #phot =Fiphot(gain=1.0, magtoflux=17., skyfit_sigma=3, skyfit_niter=4, disjoint_radius=2, apertures='5:8.0:3.0,7:10:5,8:11:6,9:12:7,10:13:7,11:14:8,12:15:9,13:16:10,14:17:11,15:18:12,20:23:17')
        phot =Fiphot(gain=1.0, magtoflux=17., skyfit_sigma=3, skyfit_niter=4, disjoint_radius=2, apertures='10:13:7,12:15:9,15:18:12,20:23:17')
        
        phot(infile, catproj, outdir=astrom_dir, outputmagnitude=True,filenum=framenb)
        #break
#exit()
cmd = "ls %s/*.fiphot > phot.in" % astrom_dir
os.system(cmd)
lc_dir = os.path.join(outpath, "LC")
lccollect_path = "lctools-collect/lctools-collect"
cmd = "%s --outdir %s" % (lccollect_path, lc_dir)
os.system(cmd)
