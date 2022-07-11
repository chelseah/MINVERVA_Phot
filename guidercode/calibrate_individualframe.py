#!/usr/bin/env python
# 
# Guiding script performing platesolving on guiding camera photos
# Should be run from MINERVAMAIN as it depends on the astrometry software to be installed
#
import os,sys
import numpy as np
import astropy.io.fits as pyfits
import glob
import time
from datetime import datetime,timedelta
import argparse
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas
import warnings
warnings.filterwarnings('ignore')

DEBUG=True

sys.path.append("/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/")
import run_astrometry_guider
astrom_dir = "/media/minervastorage/Dev/gzhou/astrom_index_5000"
cfgfile = "/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/guider_astrometry.cfg"

def calibrate(input_filename,output_filename,bias,dark,flat):
    exptime = pyfits.open(input_filename)[0].header['EXPTIME']
    fitdata = np.array(pyfits.open(input_filename)[0].data).astype(float)
    fitdata -= bias
    fitdata -= dark*(exptime/100.)
    fitdata /= flat

    fitdata += 1000

    mask = fitdata < 0
    mask += fitdata > 65000
    fitdata[mask] = np.median(fitdata.flatten())

    hdu = pyfits.open(input_filename)
    hdu[0].data = fitdata.astype(np.uint16)

    hdu.writeto(output_filename,overwrite=True)


def getlatestframe(args,bias,dark,flat):
    while True:
        #list_of_files = np.sort(glob.glob(os.path.join(args.input_dir,"*.fit")))
        #latest_file = max(list_of_files, key=os.path.getctime)
        #if datetime.now()-datetime.fromtimestamp(os.path.getmtime(latest_file)) < timedelta(days=100,minutes=1):
        while True:
            latest_file = args.input_frame
            output_file = os.path.join("T"+str(args.telescope),"caltest")
            output_file = os.path.join(output_file,os.path.basename(latest_file))

            if DEBUG:
                print("Calibrating",latest_file)
                
            calibrate(latest_file,output_file,bias,dark,flat)

            try:
                __,wcsfile,w,ra0,dec0 = run_astrometry_guider.main(output_file,os.path.join("T"+str(args.telescope),"caltest"),args.ra,args.dec,astrom_dir,cfgfile)

                ### compute current ra dec of the fiber position
                sky = w.pixel_to_world(args.fiberx,args.fibery)

                ### compute ra dec offsets to get to the fiber postion
                d_ra = float((sky.ra-args.ra*u.deg)/u.deg)
                d_dec = float((sky.dec-args.dec*u.deg)/u.deg)

                ### compute field rotation
                x0,y0 = w.world_to_pixel(SkyCoord(ra0,dec0, unit="deg"))
                x1,y1 = w.world_to_pixel(SkyCoord(ra0,dec0+0.1, unit="deg"))

                rotation = np.arctan2((x1-x0),(y1-y0))*180/np.pi
            except:
                ra0 = np.nan
                dec0 = np.nan
                d_ra = np.nan
                d_dec = np.nan
                rotation = np.nan
                
            d = {'fitsname': latest_file,
                 'file_datetime': str(datetime.fromtimestamp(os.path.getmtime(latest_file))),
                 'calc_datetime': str(datetime.now()),
                 'field_ra' : float(ra0),
                 'field_dec' : float(dec0),
                 'ra_offset' : d_ra,
                 'dec_offset' : d_dec,
                 'rotation' : rotation
            }

            if DEBUG:
                print(d)

            df = pandas.DataFrame(data=d,index=[0])
            df.to_json(os.path.join(os.path.join("T"+str(args.telescope),"caltest"),"wcssolution.json"))
            
            #sys.exit()
            #time.sleep(100)

if __name__ == "__main__":
    
    # Initiate the parser
    parser = argparse.ArgumentParser(description="Calibrate and run astrometry on Guider frames")

    # Add long and short argument
    parser.add_argument("--input_frame", "-i", help="New fits frame")
    parser.add_argument("--telescope", "-t", help="Telescope 1,3,4,5",type=int)
    parser.add_argument("--ra", help="RA in Degrees",default='NULL',type=float)
    parser.add_argument("--dec", help="DEC in Degrees",default='NULL',type=float)
    parser.add_argument("--fiberx", help="Fiber X position on detector",default='NULL',type=float)
    parser.add_argument("--fibery", help="Fiber Y position on detector",default='NULL',type=float)
    args = parser.parse_args()

    ### Check the inputs
    if args.input_frame == None:
        raise ValueError("You did not define any input frames")    

    ### get calibration frames
    cal_dir = "T"+str(args.telescope)
    bias = np.load(os.path.join(cal_dir,"bias.npy"))
    flat = np.load(os.path.join(cal_dir,"flat.npy"))
    dark = np.load(os.path.join(cal_dir,"dark.npy"))

    getlatestframe(args,bias,dark,flat)
    
    
                   
