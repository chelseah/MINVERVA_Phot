#!/usr/bin/env python
# 
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT) 
# 
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 


from __future__ import print_function
import math
import numpy as np
import scipy as sp
from scipy import linalg
import tempfile
import os

import astropy
import astropy.wcs as wcs
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord 
from astropy.time import Time
#import astroquery
#from astroquery.vizier import Vizier
import pandas as pd
from utils.dataio import readcolumn
from utils.util import ls_polynxy,poly_ordern_eval 




class Starlist(object):

    def __init__(self, infile, colid=1, colra=2, coldec=3, colmag=4,
                 colx=5, coly=6):
        self.name = infile
        self.colid = colid
        self.colra = colra
        self.coldec = coldec
        self.colmag = colmag
        self.colx = colx
        self.coly = coly
        return


class Catalog(object):
    def __init__(self, catfile, ra, dec, width, colid=1, colra=2, coldec=3,
                 colmag=4, colx=5, coly=6, ra0=None, dec0=None, method=None, maglim=14.5):
        self.name = catfile
        self.ra = float(ra)
        self.dec = float(dec)
        if ra0 is None:
            self.ra0 = self.ra
        else:
            self.ra0 = float(ra0)
        if dec0 is None:
            self.dec0 = self.dec
        else:
            self.dec0 = float(dec0)
        self.width = float(width)
        self.colid = colid
        self.colra = colra
        self.coldec = coldec
        self.colmag = colmag
        self.colx = colx
        self.coly = coly
        self.method = method
        self.maglim = float(maglim)
        print("Config Catalog: name=%s, ra=%f, dec=%f, ra0=%f, dec0=%f, width=%f, method=%s, maglim=%f" % (self.name, self.ra, self.dec, self.ra0, self.dec0, self.width, self.method, self.maglim))
    def query(self, maglim=None):
        if maglim is None:
            maglim = self.maglim
        if self.method is None:
            # query the catalog with Vizier 
            v = Vizier(columns=['UCAC4', '_RAJ2000', 'e_RAJ2000', '_DEJ2000',
                                'e_DEJ2000', 'Jmag', 'Kmag', 'Vmag', 'pmRA', 'pmDE',
                                'e_pmRA', 'e_pmDE', 'rmag', 'imag'],
                       row_limit="unlimited", column_filters={"Vmag": "<%f" % (maglim)})
            # print self.ra,self.dec
            result = v.query_region(coord.SkyCoord(ra=self.ra*u.degree,
                                                   dec=self.dec*u.degree,
                                                   frame='icrs'),
                                    width=self.width*u.degree, catalog=["UCAC4"])

            ids = result[0]['UCAC4']
            ras = result[0]['_RAJ2000'][:]
            decs = result[0]['_DEJ2000'][:]
            e_ra = result[0]['e_RAJ2000'][:]
            e_de = result[0]['e_DEJ2000'][:]
            jmag = result[0]['Jmag'][:]
            kmag = result[0]['Kmag'][:]
            vmag = result[0]['Vmag'][:]
            with open(self.name, mode='w') as fout:
                for i in range(len(ids)):
                    # ID[i]=''.join(ID[i].split('-'))
                    if jmag[i] == "--":
                        jmag[i] = np.nan
                    if kmag[i] == "--":
                        kmag[i] = np.nan
                    if vmag[i] == "--":
                        vmag[i] = np.nan
                    fout.write("%s %12.7f %12.7f %5.3f %d %d %5.3f %5.3f %5.3f\n" %
                               (ids[i], ras[i], decs[i], vmag[i], e_ra[i], e_de[i],
                                jmag[i], kmag[i], vmag[i]))
        elif self.method == 'TIC':
            # query the catalog from tic
            #TBD: reimplement a MAST TIC query
            from astroquery.mast import Catalogs
            catalog_data = Catalogs.query_region("%f %f" % (self.ra0, self.dec0), catalog="TIC", radius= self.width, Tmag=[-3, self.maglim]).to_pandas()
            #print(catalog_data.columns) 
            epoch0 = 2015.5
            epoch_now = 2021.5
            catalog_data = catalog_data[catalog_data["plx"]>0]
            c = SkyCoord(np.array(catalog_data["RA_orig"])*u.degree, np.array(catalog_data["Dec_orig"])*u.degree, pm_ra_cosdec = np.array(catalog_data["pmRA"])*np.cos(np.array(catalog_data["Dec_orig"])/180.*np.pi)*u.mas/u.yr, pm_dec=np.array(catalog_data["pmDEC"])*u.mas/u.yr, distance=(1000./np.array(catalog_data["plx"]))*u.pc, obstime=Time('2015-06-30'))
            c1 = c.apply_space_motion(new_obstime=Time('2021-06-30'))

            #newra = np.array(catalog_data["RA_orig"]) + (2021.5-2015.5) * catalog_data["pmRA"] * 1.e-3 /3600. * np.cos(np.array(catalog_data["Dec_orig"])/180.*np.pi) 
            #newdec = np.array(catalog_data["Dec_orig"]) + (2021.5-2015.5) * catalog_data["pmDEC"] * 1.e-3/3600.
            catalog_data["ra"] = c1.ra
            catalog_data["dec"] = c1.dec
            catalog_data.to_csv(self.name, columns=["ID","ra","dec","Tmag","pmRA","pmDEC","Jmag","Hmag","Kmag"], header=False, index=False)
        elif self.method == 'Gaia':
            from astroquery.mast import Catalogs
            catalog_data = Catalogs.query_region("%f %f" % (self.ra0, self.dec0), catalog="Gaia", radius= "%s deg" % (self.width), version=2).to_pandas()
            newra = np.array(catalog_data["ra"]) + (2021.5-2015.5) * catalog_data["pmRA"] * 1.e-3 /3600. * np.cos(np.array(catalog_data["dec"])/180.*np.pi) 
            newdec = np.array(catalog_data["dec"]) + (2021.5-2015.5) * catalog_data["pmDEC"] * 1.e-3/3600.
            catalog_data["ra"] = newra
            catalog_data["dec"] = newdec
            print(self.ra0, self.dec0, self.width)
            catalog_data.to_csv(self.name, columns=["source_id","ra","dec","phot_g_mean_mag","pmRA","pmdDEC","phot_bp_mean_mag","phot_gp_mean_mag","radius"], header=False, index=False)
        else:
            raise AttributeError("not implemented yet for catalog " \
                                  "query method %s" % self.method) 
        return

    




class Fistar(object):
    def __init__(self, threshold=200, sdkflag=False, fiign=False):
        super(Fistar, self).__init__()
        self.threshold = int(threshold)
        self.sdkflag = sdkflag
        self.fiign = fiign
        self.ext = '.fistar'
        
        print("Configure fistar: threshold=%f, sdkflag=%d", self.threshold, self.sdkflag) 
    def get_outfile(self, infile):
        outfile = os.path.basename(os.path.splitext(infile)[0])+self.ext
        return outfile
    def __call__(self, infile, outfile='', out_dir=''):
        if outfile == '':
            outfile = self.get_outfile(infile)
        if out_dir == '':
            out_dir = os.path.dirname(infile)
        if self.fiign:
            cmdline = 'fiign -i %s -n -o %s' % (infile, infile)
            os.system(cmdline) 
        cmdline = 'fistar "%s" -o "%s" -s flux ' \
                  '--model elliptic --flux-threshold %d --algorithm uplink ' \
                  '--iterations symmetric=2,general=1 ' \
                  '--format id,x,y,bg,amp,s,d,k,flux,s/n --comment ' \
                  % (infile, out_dir+'/'+outfile, self.threshold)
        if not self.sdkflag:
            cmdline+='--only-candidates '
        print("Excuting Fistar commands: %s", cmdline) 
        os.system(cmdline) 
        return


class Sextractor(object):
    def __init__(self, threshold=10, ext='.sexsource'):
        super(Sextractor, self).__init__()
        self.threshold = int(threshold)
        self.ext = ext
        
        print("Configure sextractor: threshold=%f", self.threshold) 
    def get_outfile(self, infile):
        outfile = os.path.basename(os.path.splitext(infile)[0])+self.ext
        return outfile
    def __call__(self, infile, outfile='', out_dir=''):
        if outfile == '':
            outfile = self.get_outfile(infile)
        if out_dir == '':
            out_dir = os.path.dirname(infile)
        cmdline = "sex %s -CATALOG_NAME %s -DETECT_THRESH %d -ANALYSIS_THRESH %d" % (infile, os.path.join(out_dir, outfile), self.threshold, self.threshold)
        print("Excuting sextractor commands: %s", cmdline) 
        os.system(cmdline)
        return



class Anet(object):
    """
    Astrometry class that uses astrometry.net to solve the frames
    """
    def __init__(self, ra=None, dec=None, q=0.01, tweak=2, order=4, radius=10, catfile=None, maxdistance=1, width=2048, height=2048, upper=0.7, lower=0.1):
        self.ra = float(ra)
        self.dec = float(dec)
        self.q = float(q)
        self.tweak = int(tweak)
        self.radius = float(radius)
        self.order = int(order)
        self.maxdistance = int(maxdistance)
        self.catfile = catfile
        self.width = width
        self.height = height
        self.upper = upper
        self.lower = lower
        # FIXME: correct the debug information
        # logger.debug("order=%d, maxdistance=%d, unitarity=%f, catfile=%s", self.order, self.maxdistance, self.unitarity, self.catfile)
         
        self.catalog = Starlist(self.catfile)
        self.cfg = None

    def get_wcsfile(self, infile): 
        wcsfile = os.path.basename(os.path.splitext(infile)[0]) + '.wcs'
        return wcsfile

    def get_transfile(self, infile): 
        transfile = os.path.basename(os.path.splitext(infile)[0]) + '.trans'
        return transfile
    
    def get_matchfile(self, infile): 
        matchfile = os.path.basename(os.path.splitext(infile)[0]) + '.match'
        return matchfile


    def get_xylsfile(self, infile):
        xylsfile = os.path.basename(os.path.splitext(infile)[0])+'.xyls'

        return xylsfile

    def polyfit(self, match_arr,transfile,xshift=0,yshift=0, totallength=0):
        i_x=match_arr[:,0]
        i_y=match_arr[:,1]
        r_x=match_arr[:,2]-xshift
        r_y=match_arr[:,3]-yshift
        if totallength ==0:
            totallength = len(i_x)
        flag=True
        reslim=0.05
        oldres=0
        #print len(i_x)
        #print np.median(i_x-r_x)
        #print np.median(i_y-r_y)
        #import matplotlib 
        #from matplotlib import pyplot as plt
        #plt.plot(i_x,i_x-r_x,'b.')
        #plt.plot(i_y,i_y-r_y,'r.')
        #plt.show()
        while flag:
            m1 = ls_polynxy(r_x, r_y, i_x, self.order)
            m2 = ls_polynxy(r_x, r_y, i_y, self.order)

            i_xe = poly_ordern_eval(r_x,r_y,m1)
            i_ye = poly_ordern_eval(r_x,r_y,m2)

            res = ((i_x-i_xe)**2.+(i_y-i_ye)**2.)**0.5
            resstd=(np.sum(res)/len(res))**0.5
            #print resstd
            if resstd<reslim or abs(resstd-oldres)/resstd<0.05:
                break
            else:
                keepindex=res<8*np.median(res)
                i_x=i_x[keepindex]
                i_y=i_y[keepindex]
                r_x=r_x[keepindex]
                r_y=r_y[keepindex]
                oldres=resstd
        #print np.median(res)
        ratio = len(i_x)*1.0/totallength
        #plt.plot(i_x,i_x-i_xe,'b.')
        #plt.plot(i_y,i_y-i_ye,'r.')
        #plt.show()
        #print m1
        #print m2
        #fout = open("/tmp/temp.txt", mode='w')
        #for i in range(len(i_x)):
        #    fout.write("%f, %f, %f, %f\n" % (i_x[i], i_y[i], r_x[i], r_y[i]))
        #fout.close()
       #print ratio
        self.output_trans(transfile, m1, m2, res, ratio)
        return
   
    def output_trans(self, transfile,m1,m2,res,ratio):
        fout = open(transfile,mode="w")
        fout.write("# Type: polynomial of order=%d (number of coefficients: %d)\n"% (self.order,len(m1)))
        fout.write("type = polynomial\n")
        fout.write("order = %d\n" % (self.order))
        fout.write("# Initial transformation of (x_img,y_img):\n")
        fout.write("offset = 0, 0\n")
        fout.write("scale = 1\n")
        fout.write("basisshift = 0, 0\n")
        fout.write("# Coefficients of the x fit:\n")
        fout.write("dxfit=")
        for i in range(len(m1)-1):
            fout.write("%14.7g," % m1[i])
        fout.write("%14.7g\n" % m1[-1])
        fout.write("# Coefficients of the y fit:\n")
        fout.write("dyfit=")
        for i in range(len(m2)-1):
            fout.write("%14.7g," % m2[i])
        fout.write("%14.7g\n" % m2[-1])
        fout.write("# Ratio: %f (percent)" % (ratio*100.))
        if len(res)<500:
            fout.write("# fit with %d stars,failed\n" % len(res))
        else:
            fout.write("# fit with %d stars\n" % len(res))
        fout.write("# Residual: %f %f\n" % ((np.sum(res)/len(res))**0.5,np.median(res)))
        fout.close()
        return
   
    def read_transfile(self, transfile):
        with open(transfile, mode='r') as fin:
            for line in fin.readlines():
                if line.startswith("order"):
                    self.order = int(line.lstrip("order = "))
                if line.startswith("dxfit"):
                    coeffs = line.lstrip("dxfit=").split(",")
                    m1 = [float(c) for c in coeffs]
                if line.startswith("dyfit"):
                    coeffs = line.lstrip("dyfit=").split(",")
                    m2 = [float(c) for c in coeffs]
        
        return [m1, m2]
    
    def kd_match(self, src, dist):
            
        import sklearn
        from sklearn.neighbors import NearestNeighbors
        nbrs=NearestNeighbors(n_neighbors=1,algorithm='auto').fit(dist)
        #distances, indices are in the shape of src
        distances,indices=nbrs.kneighbors(src)
        
        index=distances[:,0]<self.maxdistance
        match_dist=dist[indices[index]]
        
        match_src=src[index,:]
        
    
        u,ind,inv=np.unique(indices[index],return_index=True,return_inverse=True)
        
        duplicate=indices[index][ind[np.bincount(inv)>1]]
        uniq_mask=ind[np.bincount(inv)==1]
        match_distance=distances[index]
        final_dist=match_dist[uniq_mask,:,:]
        final_src=match_src[uniq_mask,:]
        
        for i in range(len(duplicate)):
            duplicate_dist=match_distance[indices[index]==duplicate[i]]
            duplicate_src=match_src[indices[index][:,0]==duplicate[i],:]
            
            keepsrc=duplicate_src[duplicate_dist==min(duplicate_dist),:]
            keepsrc=keepsrc[0].reshape([1,2])
            final_src=np.concatenate((final_src,keepsrc),axis=0)
            keepdist=dist[duplicate[i]].reshape([1,1,2])
            
            final_dist=np.concatenate((final_dist,keepdist),axis=0)
        match_arr=np.hstack([final_dist[:,0,:],final_src])    
        return match_arr




    def proj_sky_to_xy(self, xylist, inputcat='', wcsfile='', transfile='', xlim=[0, 2048], ylim=[0, 2048], out_dir=''):
        if wcsfile == '':
            wcsfile = self.get_wcsfile(xylist)
        if out_dir =='':
            out_dir = os.path.dirname(xylist)
        if inputcat == '':
            inputcat = self.catalog
        # FIXME: read input cat
        
        df = pd.read_csv(inputcat, names=["id","ra","dec", "mag", "pmra", "pmdec", "mag1", "mag2", "mag3"], index_col=False)
        ids = np.array(df["id"]).astype(str)
        ra = np.array(df.ra).astype(float)
        dec = np.array(df.dec).astype(float)
        mag = np.array(df.mag).astype(float)
        print("projecting catalogs...")
        x0, y0 = self._wcs_sky_to_xy(ra, dec, os.path.join(out_dir, wcsfile))
        if not transfile=='': 
            x, y = self._poly_xy_to_xy(x0, y0, os.path.join(out_dir, transfile))
        else:
            x = x0
            y = y0
        index = (x>xlim[0]) & (x<xlim[1]) & (y>ylim[0]) & (y<ylim[1])

        ids = ids[index]
        ra = ra[index]
        dec = dec[index]
        mag = mag[index]
        x = x[index]
        y = y[index]
        fout = open(xylist, mode="w")
        for i in range(len(ids)):
            fout.write("%s %f %f %f %f %f\n" % (ids[i], ra[i], dec[i], mag[i], x[i], y[i]))
        fout.close()
        return
    
    @staticmethod
    def _wcs_sky_to_xy(ra, dec, wcsfile):

        w = wcs.WCS(wcsfile)
        #x, y = w.all_world2pix(ra, dec, 1, maxiter=30,tolerance=1.e-3, detect_divergence=False, adaptive=True)
        x, y = w.all_world2pix(ra, dec, 1)
        #fout = open("/tmp/temp.txt", mode='a')
        #fout.write("%s\n" % wcsfile)
        #for i in range(len(x)):
        #    fout.write("%f, %f, %f, %f\n" % (ra[i], dec[i], x[i], y[i]))
        #fout.close()
       
        #x-=0.5
        #y-=0.5
        return [x, y]


    def _poly_xy_to_xy(self, x0, y0, transfile):
        m1, m2 = self.read_transfile(transfile)
        #print m1
        #print m2
        x = poly_ordern_eval(x0, y0, m1)
        y = poly_ordern_eval(x0, y0, m2)
        return [x, y]

    def __call__(self, starlist, wcsfile='', transfile='', out_dir='', matchfile='', continuesolve=False, refine=False, source='fistar', astromcfg=None):
        if wcsfile == '':
            wcsfile = self.get_wcsfile(starlist.name) 
        if transfile == '':
            transfile = self.get_transfile(starlist.name) 
        if out_dir =='':
            out_dir = os.path.dirname(starlist.name)
        if not matchfile =='':
            output_matched = True
        else:
            output_matched = False
        
        if not continuesolve:
            import tempfile
            infile = tempfile.NamedTemporaryFile(delete=False)
            if source == 'fistar':
                converting = "text2fits -H 'id X Y bg amp S D K flux SN' -f dfffffffff -n '-' '%s' '%s' " % (starlist.name, infile.name)
            elif source == 'sex':
                converting = "text2fits -H 'id flux fluxerr X Y' -f dffff -n '-' '%s' '%s'" % (starlist.name, infile.name)
            print("convert ascii to fits: %s", converting)
            os.system(converting)
            if continuesolve:
                cmdline = 'solve-field --continue "%s" -w %d -e %d  -u app --ra %f --dec %f --radius %f -p -q %.2f --tweak %d --wcs "%s" -M none -R none -B none -P none -L %f -H %f' % (infile.name, self.width, self.height, self.ra, self.dec, self.radius, self.q, self.tweak, out_dir+'/'+wcsfile, self.lower, self.upper)
            else:
                cmdline = 'solve-field --overwrite "%s" -w %d -e %d -u app --ra %f --dec %f --radius %f -p -q %.2f --tweak %d --wcs "%s" -M none -R none -B none -P none -L %f -H %f' % (infile.name, self.width, self.height, self.ra, self.dec, self.radius, self.q, self.tweak, out_dir+'/'+wcsfile, self.lower, self.upper)
            
            if not astromcfg is None:
                cmdline+=" --config %s" % astromcfg
            print("Excute astrometry.net: %s", cmdline)
            os.system(cmdline)
            os.unlink(infile.name)
        
        if refine:
            print("Project astrometry.net solution")
                
            # FIXME: read input cat
            ra = []; readcolumn(ra, self.catalog.colra, self.catalog.name); ra= np.array(ra)
            dec = []; readcolumn(dec, self.catalog.coldec, self.catalog.name); dec= np.array(dec)

            # FIXME: read input starlist 
            star_x = []; readcolumn(star_x, 2, starlist.name); star_x= np.array(star_x)
            star_y = []; readcolumn(star_y, 3, starlist.name); star_y= np.array(star_y)
            x0, y0 = self._wcs_sky_to_xy(ra, dec, os.path.join(out_dir,wcsfile))
            #print wcsfile
            index = np.isnan(star_x)
            src = np.swapaxes(np.vstack([x0, y0]), 0, 1)
            dist = np.swapaxes(np.vstack([star_x[~index], star_y[~index]]), 0, 1)
            #print src
            #print dist
            
            match_arr = self.kd_match(src,dist)
            self.polyfit(match_arr, os.path.join(out_dir, transfile), totallength=len(x0))
            
            if output_matched:
                x, y = self. _poly_xy_to_xy(match_arr[:,2], match_arr[:,3], os.path.join(out_dir, transfile))
                fmatch = open(matchfile, mode='w')
                for i in range(match_arr.shape[0]):
                    fmatch.write("%f %f %f %f\n" % (match_arr[i,0], match_arr[i,1], x[i], y[i]))
                fmatch.close()






class Fiphot(object):
    """
    photometry wrapper that interact with fiphot
    """
    def __init__(self, gain=2.1, magtoflux=21., skyfit_sigma=3, skyfit_niter=4, disjoint_radius=2, apertures='2.5:4.0:3.0'):
        super(Fiphot, self).__init__()
        self.gain = float(gain)
        self.magtoflux = float(magtoflux)
        self.skyfit_sigma = int(skyfit_sigma)
        self.skyfit_niter = int(skyfit_niter)
        self.disjoint_radius = int(disjoint_radius)
        # need to be more optimally decided
        self.apertures = apertures 
        print("fiphot configuration: gain=%.2f, magtoflux=%.1f, skyfit_sigma=%d, skyfit_niter=%d, disjoint_radius=%d, apertures=%s", self.gain, self.magtoflux, self.skyfit_sigma, self.skyfit_niter, self.disjoint_radius, self.apertures) 

    def getphotfile(self, frame):
        photfile = os.path.basename(os.path.splitext(frame)[0]) + '.fiphot'
        return photfile
        
    def __call__(self, frame, xylist, outfile='', outdir='', dry=False, positiononly=False, outputmagnitude=True, aplist=None, filenum=''):
        if outfile == '':
            outfile = self.getphotfile(frame)
        if outdir == '':
            outdir = os.path.dirname(xylist.name)
        # come up with some more general way to generate file num
        if filenum=='':
            filenum = os.path.basename(frame).split('.')[0].split('-')[1]
        if positiononly:
            # doing fiphot and then iterate on the centroid position
            cmdline = 'fiphot --input "%s" --input-list "%s" --col-id %d --col-xy %d,%d ' \
                    '--gain %.2f --mag-flux %f,10000 --apertures 4.0:5:5 ' \
                  '--sky-fit "median,sigma=%d,iterations=%d" ' \
                  '--disjoint-radius %d --format "ISXY,XY" ' \
                  '--nan-string "NaN" --aperture-mask-ignore "saturated" ' \
                  '--comment "--comment" --output - --serial %s' \
                  "| awk '{d=($5-$3)*($5-$3)+($6-$4)*($6-$4); if (0<$5 && 0<$6 && d<0.2*0.2) print $1,$5,$6;}' " \
                  '| fiphot --input "%s" --input-list - --col-id 1 --col-xy 2,3 '\
                    "--gain %.2f --mag-flux %f,10000 --apertures 4.0:5:5 " \
                  "--sky-fit 'median,sigma=%d,iterations=%d' " \
                  "--disjoint-radius %d --format 'ISXY,XY' " \
                  "--nan-string 'NaN' --aperture-mask-ignore 'saturated' " \
                  "--comment '--comment' --output - --serial %s" \
                   "| awk '{d=($5-$3)*($5-$3)+($6-$4)*($6-$4); if (0<$5 && 0<$6 && d<0.2*0.2) print $1,$5,$6;}' " \
                  '| fiphot --input "%s" --input-list - --col-id 1 --col-xy 2,3 ' \
                    "--gain %.2f --mag-flux %f,10000 --apertures 4.0:5:5 " \
                  "--sky-fit 'median,sigma=%d,iterations=%d' " \
                  "--disjoint-radius %d --format 'ISXY,XY' " \
                  "--nan-string 'NaN' --aperture-mask-ignore 'saturated' " \
                  "--comment '--comment' --output %s --serial %s" \
                  % (frame, xylist.name, xylist.colid, xylist.colx,
                     xylist.coly, self.gain, self.magtoflux, 
                     self.skyfit_sigma, self.skyfit_niter,
                     self.disjoint_radius, filenum, frame, self.gain, self.magtoflux, self.skyfit_sigma, self.skyfit_niter, self.disjoint_radius, filenum , frame, self.gain, self.magtoflux, self.skyfit_sigma, self.skyfit_niter, self.disjoint_radius, outdir+'/'+outfile, filenum)
        else:
            if outputmagnitude:
                # standard fiphot command
                if aplist is None:
                    cmdline = 'fiphot --input "%s" --input-list "%s" --col-id %d --col-xy %d,%d ' \
                          "--gain %.2f --mag-flux %f,10000 --apertures %s " \
                          "--sky-fit 'mad,sigma=%d,iterations=%d' " \
                          "--disjoint-radius %d --format 'ISXY,MmBbXYWDKws' " \
                          "--nan-string 'NaN' --aperture-mask-ignore 'saturated' " \
                          '--comment "--comment" --output "%s" --serial %s' \
                          % (frame, xylist.name, xylist.colid, xylist.colx,
                             xylist.coly, self.gain, self.magtoflux, self.apertures,
                             self.skyfit_sigma, self.skyfit_niter,
                             self.disjoint_radius, outdir+'/'+outfile, filenum)
                else:
                    cmdline = 'fiphot --input "%s" --input-list "%s" --col-id %d --col-xy %d,%d ' \
                          "--gain %.2f --mag-flux %f,10000 --apertures %s " \
                          "--sky-fit 'mad,sigma=%d,iterations=%d' " \
                          "--disjoint-radius %d --format 'ISXY,MmBbXYWDKws' " \
                          "--nan-string 'NaN' --aperture-mask-ignore 'saturated' " \
                          '--comment "--comment" --output "%s" --serial %s' \
                          % (frame, xylist.name, xylist.colid, xylist.colx,
                             xylist.coly, self.gain, self.magtoflux, aplist,
                             self.skyfit_sigma, self.skyfit_niter,
                             self.disjoint_radius, outdir+'/'+outfile, filenum)

            else:
                # output flux and do not subtract background for this simple step.
                if aplist is None:
                    cmdline = 'fiphot --input "%s" --input-list "%s" --col-id %d --col-xy %d,%d ' \
                          "--gain %.2f --apertures %s " \
                          "--sky-fit 'force=0' " \
                          "--disjoint-radius %d --format 'ISXY,FfBbXYs' " \
                          "--nan-string 'NaN' --aperture-mask-ignore 'saturated' " \
                          '--comment "--comment" --output "%s" --serial %s' \
                          % (frame, xylist.name, xylist.colid, xylist.colx,
                             xylist.coly, self.gain,  self.apertures,
                             self.disjoint_radius, os.path.join(outdir,outfile), filenum)
                else:
                    cmdline = 'fiphot --input "%s" --input-list "%s" --col-id %d --col-xy %d,%d ' \
                          "--gain %.2f --apertures %s " \
                          "--sky-fit 'force=0' " \
                          "--disjoint-radius %d --format 'ISXY,FfBbXYs' " \
                          "--nan-string 'NaN' --aperture-mask-ignore 'saturated' " \
                          '--comment "--comment" --output "%s" --serial %s' \
                          % (frame, xylist.name, xylist.colid, xylist.colx,
                             xylist.coly, self.gain,  aplist,
                             self.disjoint_radius, os.path.join(outdir,+outfile), filenum)

        print("Excute fiphot command: %s" % cmdline)
        os.system(cmdline)

