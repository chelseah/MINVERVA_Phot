#
# run_photometry_main.py controls camera and guiding loops
# Needs pywin32, astropy, requests
#

import numpy as np
import argparse
import os
import sys
from datetime import datetime
from astropy.io import fits
from astropy.time import Time
import time
import win32com.client

### import guiding functions to connect to pwi
import guiding_functions

def connect_camera():

    ### Connecting cam 1 and 2 to preserve existing ascom settings
    
    driver1 = "ASCOM.ASICamera2.Camera" 
    camera1 = win32com.client.Dispatch(driver1)
    driver2 = "ASCOM.ASICamera2_2.Camera"
    camera = win32com.client.Dispatch(driver2)

    # connect and turn on cooler
    camera.connected = True
    camera.CoolerOn = True
    camera.Offset = 30.0	
    camera.Gain = 0.
    camera.SetCCDTemperature = -5

    return camera

def take_exposure(camera,exptime,outname,exposure_settings):
    
    openshutter = True
    t = Time(datetime.utcnow(),scale='utc')
    camera.StartExposure(exptime, openshutter)
    while not camera.ImageReady: 
        time.sleep(1)	
        #print("Waiting for Image Ready")
    image = camera.ImageArray

    hdu = fits.PrimaryHDU(np.transpose(np.array(image)))
    ra,dec,alt,az,rot,focus = guiding_functions.mount_status()
    
    hdu.header['exptime'] = exptime
    hdu.header['obstime'] = t.iso
    hdu.header['jd'] = t.jd
    hdu.header['telRA'] = ra
    hdu.header['telDEC'] = dec
    hdu.header['telalt'] = alt
    hdu.header['telaz'] = az
    hdu.header['rotator'] = rot
    hdu.header['focus'] = focus
    hdu.header['CoolerOn'] = str(camera.CoolerOn)
    hdu.header['CoolerPower'] = camera.CoolerPower
    hdu.header['CCDTemp'] = camera.CCDTemperature

    try:
        hdu.header['OBJECT'] = exposure_settings['objectname']
        hdu.header['targRA'] = exposure_settings['ra']
        hdu.header['targDEC'] = exposure_settings['dec']
        hdu.header['tel'] = exposure_settings['telescope']
        hdu.header['cam'] = exposure_settings['cam']
    except Exception as e:
        print(e)
    
    hdu.writeto(outname,overwrite=True)

def close_camera(camera):
    print('turn off cooler and disconnect ')
    camera.CoolerOn = False
    camera.connected = False


def setup_tel(exposure_settings):
    print("GOING to target")

    guiding_functions.goto(exposure_settings['ra'],exposure_settings['dec'])

    while not guiding_functions.is_mount_settled() or not guiding_functions.is_rotator_settled():
        print("waiting for mount and rotator to settle")
        print("mount settled",guiding_functions.is_mount_settled())
        print("rotator settled",guiding_functions.is_rotator_settled())
        time.sleep(2)
        
    print("Mount has settled")
    

def exposure_guider_loop(camera,outname,exposure_settings):

    print("taking exposure")
    take_exposure(camera,exposure_settings["exptime"],outname,exposure_settings)

    shiftdata = guiding_functions.getshift(exposure_settings['telescope'],outname)

    if shiftdata != None:
        guider_functions.do_shift(shiftdata)

        while not guiding_functions.is_mount_settled():
            print("waiting for mount to settle")
            print("mount settled",guiding_functions.is_mount_settled())
            time.sleep(1)
        time.sleep(5)

    else:
        print("no new guiding inputs, continuing to next exposure")


def iterate(camera,outdir,exposure_settings):
    import glob
    nit = len(glob.glob(outdir))+1
    
    while nit < 5:
        
        outname = str(exposure_settings['objectname'])+"_"+str(nit)+".fits"
        outname = outdir + outname
        print("exposing ",outname)
        exposure_guider_loop(camera,outname,exposure_settings)

        nit += 1

def main(camera):

    ra = np.float(sys.argv[1])
    dec = np.float(sys.argv[2])
    telescope = 3
    objectname = "TIC124206468"
    outdir = "M:\\Data\\Photometry\\2022\\2022-07-09\\"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    exposure_settings = {
        'ra' : ra,
        'dec' : dec,
        'objectname' : objectname,
        'exptime' : 20,
        'cam' : 'ZWO294',
        'telescope' : telescope }

    outdir = outdir + "T"+str(exposure_settings['telescope'])+"\\"
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    

    setup_tel(exposure_settings)
    iterate(camera,outdir,exposure_settings)



    
if __name__ == "__main__":


    camera = connect_camera()
    
    try:
        main(camera)
    except Exception as e:
        print(e)
    except KeyboardInterrupt:
        print("Ending script gracefully")
    finally:
        close_camera(camera)

    close_camera(camera)


    
