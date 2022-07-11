import requests
import json
import os.path, time, sys
import urllib.request
import xml.etree.ElementTree as ET
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

FILE_MODIFIED_TIME_DIFF = 200 # maximum allowed age in seconds of the JSON data file 
lastfile = "" ### used to check if the current FITS file is the same as the last one


# Convert RA (deg) to H M S:
def deg2HMS( RAin ):

   if(RAin<0):
      sign = -1
      ra   = -RAin
   else:
      sign = 1
      ra   = RAin

   h = int( ra/15. )
   ra -= h*15.
   m = int( ra*4.)
   ra -= m/4.
   s = ra*240.

   if(sign == -1):
      out = '-%02d %02d %06.3f'%(h,m,s)
   else: out = '%02d %02d %06.3f'%(h,m,s)
   
   return out
   
# Convert Decl. (deg) to D M S:
def deg2DMS( Decin ):

   if(Decin<0):
      sign = -1
      dec  = -Decin
   else:
      sign = 1
      dec  = Decin

   d = int( dec )
   dec -= d
   dec *= 100.
   m = int( dec*3./5. )
   dec -= m*5./3.
   s = dec*180./5.

   if(sign == -1):
      out = '-%02d %02d %06.3f'%(d,m,s)
   else: out = '+%02d %02d %06.3f'%(d,m,s)

   return out


# Convert HH:MM:SS.SSS into Degrees :
def convHMS(ra):
   try :
      sep1 = ra.find(' ')
      hh=int(ra[0:sep1])
      sep2 = ra[sep1+1:].find(' ')
      mm=int(ra[sep1+1:sep1+sep2+1])
      ss=float(ra[sep1+sep2+2:])
   except:
      raise
   else:
      pass
   
   return(hh*15.+mm/4.+ss/240.)

# Convert +DD:MM:SS.SSS into Degrees :
def convDMS(dec):

   Csign=dec[0]
   if Csign=='-':
      sign=-1.
      off = 1
   elif Csign=='+':
      sign= 1.
      off = 1
   else:
      sign= 1.
      off = 0

   try :
      sep1 = dec.find(' ')
      deg=int(dec[off:sep1])
      sep2 = dec[sep1+1:].find(' ')
      arcmin=int(dec[sep1+1:sep1+sep2+1])
      arcsec=float(dec[sep1+sep2+2:])
   except:
      raise
   else:
      pass

   return(sign*(deg+(arcmin*5./3.+arcsec*5./180.)/100.))



def getshift(tel,lastfile):

    data = None
    
    wcsfile = 'M:\\Dev\\gzhou\\temp\\wcssolution_'+str(int(tel))+'.json'
    
    mod_time = os.path.getmtime(wcsfile) 
    secs_elapsed = time.time() - mod_time

    if ( round(secs_elapsed) < FILE_MODIFIED_TIME_DIFF) :

        print("Reading new JSON file...")
        
        with open(p) as f:
            data = json.load(f)
            if lastfile == data['fitsname']['0']:
                data = None
    return data


def do_shift(data):

    ra_offset_arcsec = float(data["ra_offset"]["0"]) * 3600
    dec_offset_arcsec = float(data["dec_offset"]["0"]) * 3600
    rot_offset = float(data["theta_offset"]["0"])

    ### scale scope settle time as multiple of exposure time
    SCOPE_SETTLE_TIME = float(data["exptime"]["0"])

    if (ra_offset_arcsec**2+dec_offset_arcsec**2) < 1.**2:
        SCALE_FACTOR = 0
    if (ra_offset_arcsec**2+dec_offset_arcsec**2) > 1.**2 and (ra_offset_arcsec**2+dec_offset_arcsec**2) < 5.**2:
        SCALE_FACTOR = 0.5
    if (ra_offset_arcsec**2+dec_offset_arcsec**2) > 5.**2:
        SCALE_FACTOR = 1

    ra_offset_arcsec *= SCALE_FACTOR
    dec_offset_arcsec *= SCALE_FACTOR

    if abs(rot_offset) > 5:
        rot_offset = 0

    rot_offset *= 0.5

    # Call PWI's Http interface to pass deltas to mount
    print("Updating mount position with ra_offset=%s  and  dec_offset=%s" % (ra_offset_arcsec,dec_offset_arcsec))
    url = "http://localhost:8080/?device=mount&cmd=move&incrementra=%s&incrementdec=%s" % (ra_offset_arcsec, dec_offset_arcsec)
    resp = requests.get(url)

    print("updating rotator position with rot offset=%s" % (rot_offset))
    url = "http://localhost:8080/?device=rotator2&cmd=move&increment=%s" % (rot_offset)
    resp = requests.get(url)

    lastfile = data['fitsname']['0']


    return lastfile,ra_offset_arcsec**2+dec_offset_arcsec**2,rot_offset



def goto(ra,dec):

    rastr = deg2HMS(ra)
    decstr = deg2DMS(dec)

    # Call PWI's Http interface to pass deltas to mount
    print("GOTO ra=%s  and  dec=%s" % (rastr,decstr))
    url = "http://localhost:8080/?device=mount&cmd=move&ra2000=%s&dec2000=%s" % (rastr,decstr)

    print(url)
    resp = requests.get(url)

    print("Starting rotator position at 150deg")
    url = "http://localhost:8080/?device=rotator2&cmd=move&position=%s" % (150)
    resp = requests.get(url)
    




def xmltodic(xmldata):
    data = ET.ElementTree(ET.fromstring(xmldata))
    root = data.getroot()

    outdic = {}
    for i in range(len(root)):
        name1 = root[i].tag
        for sts in root[i]:
            name2 = sts.tag
            value = sts.text

            outdic['%s_%s' % (name1,name2)] = value

    return outdic

def mount_status():
    url = "http://localhost:8080/?cmd=getsystem"    
    resp = urllib.request.urlopen(url).read()
    mtstatus = xmltodic(resp)
    
    ra = mtstatus['mount_ra_2000']
    dec = mtstatus['mount_dec_2000']
    alt = float(mtstatus['mount_alt_radian'])*180/np.pi
    az = float(mtstatus['mount_azm_radian'])*180/np.pi

    try:
        rot = mtstatus['rotator2_position']
        focus = mtstatus['focuser2_position']
    except Exception as e:
        print(e)
        rot = -99
        focus = -99
    
    return ra,dec,alt,az,rot,focus


def is_mount_settled():
    url = "http://localhost:8080/?cmd=getsystem"    
    resp = urllib.request.urlopen(url).read()
    mtstatus = xmltodic(resp)

    ra = convHMS(mtstatus['mount_ra'])
    dec = convDMS(mtstatus['mount_dec'])

    reqra = convHMS(mtstatus['mount_ra_target'])
    reqdec = convDMS(mtstatus['mount_dec_target'])

    dist = np.sqrt(((ra-reqra)*np.cos(dec*np.pi/180))**2 + (dec-reqdec)**2)

    print("Distance from target ",dist,ra,reqra,dec,reqdec)

    return dist < 0.02

    
def is_rotator_settled():
    url = "http://localhost:8080/?cmd=getsystem"    
    resp = urllib.request.urlopen(url).read()
    mtstatus = xmltodic(resp)
    
    return eval(mtstatus['rotator2_goto_complete'])
