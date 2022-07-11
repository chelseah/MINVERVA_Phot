import requests
import json
import os.path, time, sys


# path to platesolving results JSON data file 
p = 'M:\\Dev\\gzhou\\temp\\wcssolution_4.json'
#p = '.\wcssolution.json'

FILE_MODIFIED_TIME_DIFF = 200 # maximum allowed age in seconds of the JSON data file 
SLEEP_TIME = 5 # time in secs to wait between each data file readings
#SCOPE_SETTLE_TIME = 5 # time in secs to wait for scope to settle after position change
SCALE_FACTOR = 1 ### Scale factor for the offset correction -- do not apply the entire prescribed correction
TIMEOUT_CUTOFF = 600 ### maximum timeout before killing loop

lastfile = "" ### used to check if the current FITS file is the same as the last one



cumulative_zeroinput_time = 0 ### counts up everytime it gets a dRa=0 dDec=0 command

# begin endless loop
while (True):
    
    try:
        #print("checking for new JSON file...")

        mod_time = os.path.getmtime(p) # time.ctime()
        secs_elapsed = time.time() - mod_time
        #hours, rest = divmod(secs_elapsed, 3600)
        #minutes, seconds = divmod(rest, 60)

        if ( round(secs_elapsed) > FILE_MODIFIED_TIME_DIFF) :
            #raise ValueError("JSON file too old")
            print("JSON file is %s seconds old !!" % (round(secs_elapsed)))
            print("** Waiting for newer file...")
            time.sleep(SLEEP_TIME)
            continue

        print("Reading new JSON file...")

        with open(p) as f:
            data = json.load(f)
    
        ra_offset_arcsec = float(data["ra_offset"]["0"]) * 3600
        dec_offset_arcsec = float(data["dec_offset"]["0"]) * 3600
        rot_offset = float(data["theta_offset"]["0"])
		
        ### scale scope settle time as multiple of exposure time
        SCOPE_SETTLE_TIME = float(data["exptime"]["0"])
        if ra_offset_arcsec**2+dec_offset_arcsec**2 < 0.7**2:
            SCOPE_SETTLE_TIME = 1
            SCALE_FACTOR = 0
        if ra_offset_arcsec**2+dec_offset_arcsec**2 > 0.7**2 and ra_offset_arcsec**2+dec_offset_arcsec**2 < 3**2:
            SCOPE_SETTLE_TIME *= 3
            SCALE_FACTOR = 0.5
        if ra_offset_arcsec**2+dec_offset_arcsec**2 > 3**2:
            SCOPE_SETTLE_TIME *= 5
            SCALE_FACTOR = 1

        ra_offset_arcsec *= SCALE_FACTOR
        dec_offset_arcsec *= SCALE_FACTOR
		
		
        #if abs(rot_offset) < 0.2:
        #    rot_offset = 0
        #else:
        #    rot_offset *= 0.7
			
        #if SCOPE_SETTLE_TIME < 10: ### set minimum scope settle time
        #    SCOPE_SETTLE_TIME = 10
		
        #if abs(rot_offset) > 5:
        #    rot_offset = 0

        ### Check if we have already applied the corrections from this FITS frame
        if lastfile == data['fitsname']['0']:
            ra_offset_arcsec = 0
            dec_offset_arcsec = 0
            rot_offset = 0
            SCOPE_SETTLE_TIME = 1
            print("Already applied current solution")

        if ra_offset_arcsec == 0 and dec_offset_arcsec == 0:
            cumulative_zeroinput_time += SLEEP_TIME

        else:
            cumulative_zeroinput_time = 0

        if cumulative_zeroinput_time > TIMEOUT_CUTOFF:
            print("Guider loop timeout. Exceeded",TIMEOUT_CUTOFF)
            sys.exit()


        # Call PWI's Http interface to pass deltas to mount
        print("Updating mount position with ra_offset=%s  and  dec_offset=%s" % (ra_offset_arcsec,dec_offset_arcsec))
        url = "http://localhost:8080/?device=mount&cmd=move&incrementra=%s&incrementdec=%s" % (ra_offset_arcsec, dec_offset_arcsec)
        resp = requests.get(url)

        print("updating rotator position with rot offset=%s" % (rot_offset))
        url = "http://localhost:8080/?device=rotator2&cmd=move&increment=%s" % (rot_offset)
        resp = requests.get(url)

        print("Waiting for scope to settle... "+str(SCOPE_SETTLE_TIME)+"s")

        time.sleep(SCOPE_SETTLE_TIME)
        print("Done")
        lastfile = data['fitsname']['0']

    #except ValueError as e:
    #    print(e.message)
    #    pass
    except Exception as e:
        print(e)
        pass

    time.sleep(SLEEP_TIME)