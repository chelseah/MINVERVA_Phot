#!/usr/bin/env python3.6
# 
# Guiding script performing platesolving on ZWO294 Science cameras
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

# datetime object containing current date and time
now = datetime.now()
dt_string = now.strftime("%d%m%Y-%H:%M:%S")

DEBUG=True

sys.path.append("/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/")
import run_astrometry_guider as run_astrometry_photometry
cfg_file = "/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/guider_astrometry_fine.cfg"
cfg_file_fine = "/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/guider_astrometry_fine.cfg"
tempdir = "/media/minervastorage/Dev/gzhou/temp"



def run_frame(args,theta_old=np.nan):

    ### Run astrometry on a single frame
    
    latest_file = args.input_frame
    output_file = os.path.join(tempdir,os.path.basename(args.input_frame))
    time.sleep(5)
    os.system("cp "+latest_file+" "+output_file)
    time.sleep(5)

    fitsheader = pyfits.open(output_file)
    ra = float(fitsheader[0]['targRA'])
    dec = float(fitsheader[0]['targDEC'])
    tel = str(fitsheader[0]['tel'])
    cam = str(fitsheader[0]['cam'])

    
    exptime = 1
    
    try:                

        #fitshape = np.shape(pyfits.open(output_file)[0].data)

        if cam == 'QHY600':
            cfg_file = "/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/cfg4100.cfg"

            __,wcsfile,w,ra0,dec0,theta = run_astrometry_photometry.main(output_file,tempdir,ra,dec,tempdir,cfg_file,cam)

            if ra0 > 1000 or dec0 > 1000:
                cfg_file = "/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/cfg4100.cfg"
                __,wcsfile,w,ra0,dec0,theta = run_astrometry_photometry.main(output_file,tempdir,ra,dec,tempdir,cfg_file_fine,cam)

        if cam == 'andor':
            cfg_file = "/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/cfg4100.cfg"

            fits = pyfits.open(args.input_frame)
            exptime = float(fits[0].header['EXPTIME'])

            __,wcsfile,w,ra0,dec0,theta = run_astrometry_photometry.main(output_file,tempdir,ra,dec,tempdir,cfg_file,cam)

        if cam == 'ZWO294':
            cfg_file = "/media/minervastorage/Dev/gzhou/MAQLP/QLP-MINERVA/bin/cfg5000.cfg"

            __,wcsfile,w,ra0,dec0,theta = run_astrometry_photometry.main(output_file,tempdir,ra,dec,tempdir,cfg_file,'ZWO295')

            fits = pyfits.open(args.input_frame)
            exptime = float(fits[0].header['EXPTIME'])

                        
        ### compute ra dec offsets to get to the fiber postion

        if DEBUG:
            print(ra0,type(ra0))
            print(dec0,type(dec0))
            print('theta',theta)


        d_ra = -1*float(ra0-ra)
        d_dec = -1*float(dec0-dec)

        if theta_old == theta_old and theta == theta:
            dtheta = theta - theta_old
        else:
            dtheta = 0

        if d_ra > 1 or d_dec > 1:
            d_ra = 0
            d_dec = 0
            dtheta = 0

    except FileNotFoundError:
        ra0 = np.nan
        dec0 = np.nan
        d_ra = 0
        d_dec = 0
        dtheta = 0


    if np.isnan(d_ra) or np.isnan(d_dec) or np.abs(d_ra) > 2 or np.abs(d_dec) > 2:
        d_ra = 0
        d_dec = 0
        dtheta = 0

    d = {'fitsname': args.input_frame,
         'file_datetime': str(datetime.fromtimestamp(os.path.getmtime(args.input_frame))),
         'calc_datetime': str(datetime.now()),
         'field_ra' : float(ra0),
         'field_dec' : float(dec0),
         'ra_offset' : d_ra,
         'dec_offset' : d_dec,
         'theta_offset' : dtheta,
         'exptime' : exptime,
    }

    if DEBUG:
        print(d)

    df = pandas.DataFrame(data=d,index=[0])
    df.to_json(os.path.join(tempdir,"wcssolution_"+str(int(args.tel))+".json"))
    return theta,d_ra,d_dec,dtheta


def runfolder(args):
    ### run astrometry on the latest image in a folder 
    
    os.system("rm "+tempdir+"/*.fits*")
    #logpath = os.path.join(tempdir,"TEL"+str(int(args.tel))+"_"+dt_string+".log")


    log = []

    theta_old = np.nan

    lastfile = ""
    while True:
        list_of_files = list(set(glob.glob(os.path.join(args.input_folder,"*.fit*")))-set(glob.glob(os.path.join(args.input_folder,"*FLAT*")))-set(glob.glob(os.path.join(args.input_folder,"*BIAS*")))-set(glob.glob(os.path.join(args.input_folder,"*DARK*"))))

        #print(list_of_files)
        if len(list_of_files) > 0:

            latest_file = max(list_of_files, key=os.path.getmtime)
            #if DEBUG:
            #    print("running ",latest_file)
            #    #sys.exit()
            if latest_file != lastfile and datetime.now()-datetime.fromtimestamp(os.path.getmtime(latest_file)) < timedelta(days=300,minutes=1):
                if DEBUG:
                    print("running",latest_file)
                
                args.input_frame = latest_file
                theta_new,d_ra,d_dec,d_theta = run_frame(args,theta_old=theta_old)

                if theta_old != theta_old:
                    theta_old = theta_new

                ### save log
                log.append([d_ra,d_dec,d_theta,theta_old])
                lastfile = latest_file

                
                #np.savetxt(logpath,np.array(log),fmt="%.10f")
        else:
            time.sleep(10)
                
            
if __name__ == "__main__":
    
    # Initiate the parser
    parser = argparse.ArgumentParser(description="Calibrate and run astrometry on Guider frames")

    # Add long and short argument
    parser.add_argument("--input_frame", "-f", help="Run astrometry on input frame")
    parser.add_argument("--input_folder", "-f", help="Run astrometry continuously on the latest image in a folder")
    args = parser.parse_args()

    ### Check the inputs
    if args.input_frame == None and args.input_folder == None:
        raise ValueError("You did not define any input frames or folders")    

    if args.input_folder != None: 
        runfolder(args)
    
    if args.input_frame != None and args.input_folder == None: 
        run_frame(args)
    
                   
