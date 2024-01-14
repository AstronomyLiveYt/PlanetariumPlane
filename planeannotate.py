import numpy as np
import datetime
from datetime import timedelta
import ephem
import math
import time
import base64
import os
import cv2
import pandas as pd
import sys
import pytz
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io.fits import getheader
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import FK5
from PIL import Image

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))


line1 = "ISS (ZARYA)"
line2 = "1 25544U 98067A   21118.05230159  .00003481  00000-0  71408-4 0  9998"
line3 = "2 25544  51.6436 230.7008 0002607 287.5877 248.8646 15.48957981280821"
iss = ephem.readtle(line1, line2, line3)

def equatorial_to_horizon(dec, ra, lat, lst):
    lst = math.radians(lst)
    ra = math.radians(ra)
    dec = math.radians(dec)
    lat = math.radians(lat)
    hour = lst - ra
    if math.degrees(hour) < 0:
        hour = math.radians(math.degrees(hour) + 360)
    elif math.degrees(hour) > 360:
        hour = math.radians(math.degrees(hour) - 360)
    alt = math.asin(math.sin(dec)*math.sin(lat)+math.cos(dec)*math.cos(lat)*math.cos(hour))
    az = math.acos((math.sin(dec)-math.sin(lat)*math.sin(alt))/(math.cos(lat)*math.cos(alt)))
    if math.sin(hour)>0:
        az = (360 - (az * 180/math.pi))*(math.pi/180)
    az = math.radians(math.degrees(az))
    if math.degrees(az) > 360:
        az = math.radians(math.degrees(az) - 360)
    return(alt, az)

font = cv2.FONT_HERSHEY_COMPLEX_SMALL
fontScale = 0.8
blankimage = np.zeros(((1080),(1080),3), np.uint8)
dnowstart = datetime.datetime(2021, 4, 30, 0, 48, 48).replace(tzinfo=datetime.timezone.utc)
print(dnowstart)

Lat = 0.0
Lon = 0.0
observatory = ephem.Observer()

trajdf = pd.read_csv("GRUtoLISflightdatav7.csv")

realplaneras = []
realplanedecs = []
coordsequence = []
expectingsequence = True

fps = 29.97
fourcc = cv2.VideoWriter_fourcc(*str('mp4v'))
out = cv2.VideoWriter('annotatedplane.mp4',fourcc, fps, ((1920),(650)))
delta = 0
for t in range(0,1689):
    delta=(t*6.764640909090909)
    framenumber = int(1000+t)
    framename = str('GRUtoLIS_'+str(framenumber)+'_output.jpg')
    annotateframename = str('GRUtoLIS_'+str(framenumber)+'_output-ngc.png')
    fitsframename = str('GRUtoLIS_'+str(framenumber)+'_output.new')
    try:
        hdu = fits.open(fitsframename)[0]
        hdu.header['NAXIS'] = 2
        wcs = WCS(hdu)
    except:
        eightbit = cv2.imread(framename, cv2.IMREAD_COLOR)
        out.write(eightbit)
        cv2.imshow('Output',eightbit)
        cv2.waitKey(1)
        continue
    dnow = dnowstart + datetime.timedelta(seconds=delta)
    #Now find corresponding time in flight trajectory data
    FoundCorrectRow = False
    for index, row in trajdf.iterrows():
        trajtime = datetime.datetime.strptime(row['Time (UT)'], "%m/%d/%Y %H:%M:%S").replace(tzinfo=datetime.timezone.utc)
        if dnow <= trajtime and FoundCorrectRow is False:
            #This row is some time in the future so we need to interpolate
            FoundCorrectRow = True
            Lat2 = float(row['Latitude'])
            Lon2 = float(row['Longitude'])
            Course2 = float(row['Course'])
            LatDelta = Lat2 - LastLat
            LonDelta = Lon2 - LastLon
            CourseDelta = Course2 - LastCourse
            RowDiff = (trajtime - LastTime).total_seconds()
            currentdt = (dnow - LastTime).total_seconds()
            percenttonext = (currentdt/RowDiff)
            CurrentLat = (LatDelta*percenttonext)+LastLat
            CurrentLon = (LonDelta*percenttonext)+LastLon
            CurrentCourse = (CourseDelta*percenttonext)+LastCourse
        LastTime = trajtime
        LastLat = float(row['Latitude'])
        LastLon = float(row['Longitude'])
        LastCourse = float(row['Course'])
    observatory.lat = str(CurrentLat)
    observatory.lon = str(CurrentLon)
    observatory.date = dnow
    alt = math.radians(5)
    az = math.radians(CurrentCourse)
    zenith = 90-math.degrees(alt)
    ra, dec = observatory.radec_of(az, alt)
    ra = ephem.degrees(ra)
    dec = ephem.degrees(dec)
    
    coordinates = [str(str(ra) + ' ' + str(dec))]
    coordinates2 = SkyCoord(coordinates, frame=FK5, unit="deg")
    pixels = skycoord_to_pixel(coordinates2,wcs,origin=0,mode='wcs')
    cx, cy = pixels
    cx = round(cx[0])
    cy = round(cy[0])
    center_coordinates = (int(cx), int(cy))
    color2 = (0,0,255)
    thickness = 5
    radius = 20
    try:
        eightbit = cv2.imread(annotateframename, cv2.IMREAD_COLOR)
    except:
        continue
    eightbit = cv2.circle(eightbit, center_coordinates, radius, color2, thickness)
    #Muscida Zenith Distance   
    lst = math.degrees(observatory.sidereal_time())
    lat = CurrentLat
    ra = 127.566125
    dec = 60.71816822222222
    alt, az = equatorial_to_horizon(dec, ra, lat, lst)
    alt = math.degrees(alt)
    az = math.degrees(az)
    zenithdist = round(90-alt,2)
    coordinates = [str(str(ra) + ' ' + str(dec))]
    coordinates2 = SkyCoord(coordinates, frame=FK5, unit="deg")
    pixels = skycoord_to_pixel(coordinates2,wcs,origin=0,mode='wcs')
    cx, cy = pixels
    cx = round(cx[0])
    cy = round(cy[0])
    center_coordinates = (int(cx), int(cy))
    color2 = (0,255,0)
    thickness = 2
    radius = 10
    org=(50,80)
    eightbit = cv2.circle(eightbit, center_coordinates, radius, color2, thickness)
    thickness = 1
    eightbit = cv2.putText(eightbit, str('Muscida Zenith Distance: '+str(zenithdist)), org, font, fontScale, color2, thickness, cv2.LINE_AA) 
    #Alderamin Zenith Distance
    ra = 319.6448847083333
    dec = 62.58557447222222
    alt, az = equatorial_to_horizon(dec, ra, lat, lst)
    alt = math.degrees(alt)
    az = math.degrees(az)
    zenithdist = round(90-alt,2)
    coordinates = [str(str(ra) + ' ' + str(dec))]
    coordinates2 = SkyCoord(coordinates, frame=FK5, unit="deg")
    pixels = skycoord_to_pixel(coordinates2,wcs,origin=0,mode='wcs')
    cx, cy = pixels
    cx = round(cx[0])
    cy = round(cy[0])
    center_coordinates = (int(cx), int(cy))
    color2 = (255,100,100)
    thickness = 2
    radius = 10
    org=(50,100)
    eightbit = cv2.circle(eightbit, center_coordinates, radius, color2, thickness)
    thickness = 1
    eightbit = cv2.putText(eightbit, str('Alderamin Zenith Distance: '+str(zenithdist)), org, font, fontScale, color2, thickness, cv2.LINE_AA) 
    #Write lat and lon
    fontScale = 1
    thickness = 1
    color2 = (255,255,255)

    latprint = round(CurrentLat,2)
    lonprint = round(CurrentLon,2)
    datestr = str(str(observatory.date)+ ' UT')
    latstr = str('Lat: '+str(latprint)+' N')
    lonstr = str('Lon: '+str(lonprint)+' E')
    org=(50,20)
    eightbit = cv2.putText(eightbit, str(datestr), org, font,  
                   fontScale, color2, thickness, cv2.LINE_AA) 
    org=(50,40)
    eightbit = cv2.putText(eightbit, str(latstr), org, font,  
                   fontScale, color2, thickness, cv2.LINE_AA) 
    org=(50,60)               
    eightbit = cv2.putText(eightbit, str(lonstr), org, font,  
                   fontScale, color2, thickness, cv2.LINE_AA) 
    #Label ISS
    try:
        iss.compute(observatory)
        issra = ephem.degrees(iss.ra)
        issdec = ephem.degrees(iss.dec)
        coordinates = [str(str(issra) + ' ' + str(issdec))]
        coordinates2 = SkyCoord(coordinates, frame=FK5, unit="deg")
        pixels = skycoord_to_pixel(coordinates2,wcs,origin=0,mode='wcs')
        cx, cy = pixels
        cx = round(cx[0])
        cy = round(cy[0])
        org = (int(cx), int(cy))
        color2 = (255,255,255)
        thickness = 1
        eightbit = cv2.putText(eightbit, str('ISS'), org, font, fontScale, color2, thickness, cv2.LINE_AA)
    except:
        pass
                   
    out.write(eightbit)
    cv2.imshow('Output',eightbit)
    cv2.waitKey(1)
out.release()
