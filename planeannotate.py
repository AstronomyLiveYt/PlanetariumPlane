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

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
constellationsfile = os.path.join(THIS_FOLDER, 'constellationlines.csv')
starfile = os.path.join(THIS_FOLDER, 'bsc5.dat')

def equatorial_to_horizon(dec, ra, lat, lst):
    hour = lst - ra
    if math.degrees(hour) < 0:
        hour = math.radians(math.degrees(hour) + 360)
    elif math.degrees(hour) > 360:
        hour = math.radians(math.degrees(hour) - 360)
    alt = math.asin(math.sin(dec)*math.sin(lat)+math.cos(dec)*math.cos(lat)*math.cos(hour))
    az = math.acos((math.sin(dec)-math.sin(lat)*math.sin(alt))/(math.cos(lat)*math.cos(alt)))
    if math.sin(hour)>0:
        az = (360 - (az * 180/math.pi))*(math.pi/180)
    az = math.radians(math.degrees(az)-90)
    if math.degrees(az) > 360:
        az = math.radians(math.degrees(az) - 360)
    return(alt, az)

O = (255, 180, 157)
B = (255, 188, 167)
A = (255, 209, 192)
F = (255, 232, 228)
G = (236, 245, 255)
K = (174, 215, 255)
M = (123, 187, 255)
font = cv2.FONT_HERSHEY_COMPLEX_SMALL
blankimage = np.zeros(((1080),(1080),3), np.uint8)
dnowstart = datetime.datetime(2016, 8, 8, 22, 00, 11).replace(tzinfo=datetime.timezone.utc)
print(dnowstart)
df = pd.read_csv(constellationsfile)

Lat = 47.4849
Lon = 8.5294
observatory = ephem.Observer()

trajdf = pd.read_csv("swr92flightdatafixedv2.csv")

astrometryfile = open('astrometrylog.txt', 'r')
realplaneras = []
realplanedecs = []
coordsequence = []
expectingsequence = True
#Pre-load astrometry coordinates from real time lapse pictures
with open('astrometrylog.txt', encoding="utf8") as f:
    lines = [line.rstrip('\n') for line in f]
    for idx, line in enumerate(lines):
        if 'solved: writing to file ./flightlapse_' in line:
            numbertuple = line.split('solved: writing to file ./flightlapse_')[1].split('_output.solved')[0]
            coordsequence.append(int(numbertuple))
            expectingsequence = False
        if 'Field center: (RA,Dec) = (' in line and expectingsequence is False:
            expectingsequence = True
            coordtuple = line.split('(RA,Dec) = (')[1].split(') deg.')[0]
            realra = math.radians(float(coordtuple.split(',')[0]))
            realdec = math.radians(float(coordtuple.split(',')[1]))
            realplaneras.append(realra)
            realplanedecs.append(realdec)

data = {'Sequence':coordsequence,'Ra':realplaneras,'Dec':realplanedecs}
astrometrydf = pd.DataFrame(data)
astrometrydf = astrometrydf.sort_values(by=['Sequence'])

fps = 25
fourcc = cv2.VideoWriter_fourcc(*str('mp4v'))
out = cv2.VideoWriter('planetariumplane.mp4',fourcc, fps, ((1080),(1080)))
delta = 0
with open(starfile) as f:
    lines = [line.rstrip('\n') for line in f]
    #for timedelta in range(0,1):
    #while True:
    for idx, realplanerow in astrometrydf.iterrows():
        delta+=5
        blankimage = np.zeros(((1080),(1080),3), np.uint8)
        dnow = dnowstart + datetime.timedelta(seconds=delta)
        #Now find corresponding time in flight trajectory data
        maxtimedelta = timedelta(days = 7)
        for index, row in trajdf.iterrows():
            trajtime = datetime.datetime.strptime(row['Time (EDT)'], "%m/%d/%Y %H:%M").replace(tzinfo=datetime.timezone.utc)
            trajtime = trajtime + datetime.timedelta(hours=4)
            currentdt = abs(dnow - trajtime)
            if currentdt < maxtimedelta:
                maxtimedelta = currentdt
                Lat = float(row['Latitude'])
                Lon = float(row['Longitude'])
                Course = float(row['Course'])
                dbest = trajtime
        observatory.lat = str(Lat)
        observatory.lon = str(Lon)
        observatory.date = dnow
        alt = math.radians(5)
        az = math.radians(Course-90)
        zenith = 90-math.degrees(alt)
        x = (zenith*math.cos(az))*(540/90)+540
        y = (zenith*math.sin(az))*(540/90)+540
        cv2.circle(blankimage,(int(x),int(y)), int(20), (0,255,0), 1)
        realra = realplanerow['Ra']        
        realdec = realplanerow['Dec']
        alt,az = equatorial_to_horizon(realdec, realra, observatory.lat, observatory.sidereal_time())
        #print(dbest, dnow, math.degrees(alt),math.degrees(az))
        if math.degrees(alt) > 00:
            zenith = 90-math.degrees(alt)
            x = (zenith*math.cos(az))*(540/90)+540
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(20), (0,0,255), 1)
        #dnow = datetime.datetime.utcnow()
        observatory.date = dnow
        #print(dnow)
        #Do constellation lines
        trueindex = -1
        for index, row in df.iterrows():
            trueindex +=1
            rah = row['ra']
            ra = math.radians(float(rah)*15)
            dec = math.radians(row['dec'])
            alt, az = equatorial_to_horizon(dec, ra, observatory.lat, observatory.sidereal_time())
            if trueindex > 0:
                try:
                    az1 = az
                    zenith = 90-math.degrees(alt)
                    x1 = int((zenith*math.cos(az1))*(540/90)+540)
                    y1 = int((zenith*math.sin(az1))*(540/90)+540)
                    
                    az2 = azlast
                    zenith = 90-math.degrees(altlast)
                    x2 = int((zenith*math.cos(az2))*(540/90)+540)
                    y2 = int((zenith*math.sin(az2))*(540/90)+540)
                    if alt > 0 and altlast > 0:
                        cv2.line(blankimage,(x1,y1),(x2,y2),(255,255,255),1)
                except:
                    trueindex = -1
                    pass
            
            altlast = alt
            azlast = az
        #Load up the stars
        for idx, line in enumerate(lines):
            rahour = line[75:77]
            ramin = line[77:79]
            rasec = line[79:83]
            decdeg = line[83:86]
            decmin = line[86:88]
            decsec = line[88:90]
            spectral = line[129:130]
            mag = line[102:107]
            try:
                pixsize = int(5.0/(2+float(mag)))+1 
                ra = str(rahour+':'+ramin+':'+rasec)
                dec = str(decdeg+':'+decmin+':'+decsec)
                star = ephem.Equatorial(ra, dec, epoch='2000')
                ra = star.ra
                dec = star.dec
                if spectral == 'O':
                    starcolor = O
                if spectral == 'B':
                    starcolor = B
                if spectral == 'A':
                    starcolor = A
                if spectral == 'F':
                    starcolor = F
                if spectral == 'G':
                    starcolor = G
                if spectral == 'K':
                    starcolor = K
                if spectral == 'M':
                    starcolor = M
                alt,az = equatorial_to_horizon(dec, ra, observatory.lat, observatory.sidereal_time())
                if math.degrees(alt) > 00 and float(mag) < 5:
                    zenith = 90-math.degrees(alt)
                    x = (zenith*math.cos(az))*(540/90)+540
                    y = (zenith*math.sin(az))*(540/90)+540
                    cv2.circle(blankimage,(int(x),int(y)), int(pixsize), starcolor, -1)
                    #print(x, y)
                    #print(math.degrees(alt), math.degrees(az))
                    #print(line)
                    #print(rahour, ramin, rasec, decdeg, decmin, decsec, spectral, mag, int(pixsize), math.degrees(alt), math.degrees(az), int(x), int(y))
            except:
                pass
        #Load up the moon, sun, and planets
        blankimage = cv2.flip(blankimage, 1)
        celestial = ephem.Moon(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(15), (200,200,200), -1)
            cv2.putText(blankimage,'Moon',(int(x+10),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        celestial = ephem.Sun(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(15), (0,255,255), -1)
            cv2.putText(blankimage,'Sun',(int(x+10),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        celestial = ephem.Mercury(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(5), (255,100,100), -1)
            cv2.putText(blankimage,'Mercury',(int(x+5),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        celestial = ephem.Venus(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(5), (255,255,255), -1)
            cv2.putText(blankimage,'Venus',(int(x+5),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        celestial = ephem.Mars(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(5), (50,50,255), -1)
            cv2.putText(blankimage,'Mars',(int(x+5),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        celestial = ephem.Jupiter(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(5), (150,150,255), -1)
            cv2.putText(blankimage,'Jupiter',(int(x+5),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        celestial = ephem.Saturn(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(5), (100,100,255), -1)
            cv2.putText(blankimage,'Saturn',(int(x+5),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        celestial = ephem.Uranus(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(5), (255,50,50), -1)
            cv2.putText(blankimage,'Uranus',(int(x+5),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        celestial = ephem.Neptune(observatory)
        alt = celestial.alt
        if math.degrees(alt) > 0:
            az = math.radians(math.degrees(celestial.az)-90)
            zenith = 90-math.degrees(alt)
            x = 540-(zenith*math.cos(az))*(540/90)
            y = (zenith*math.sin(az))*(540/90)+540
            cv2.circle(blankimage,(int(x),int(y)), int(5), (255,0,0), -1)
            cv2.putText(blankimage,'Neptune',(int(x+5),int(y)), font, 1,(255,255,255),1,cv2.LINE_AA)
        cv2.circle(blankimage,(int(540),int(540)), int(540), (255,255,255), 1)
        tlefile = ""
        datetimestringprint = str(dnow.strftime('%m-%d-%Y %H:%M:%S'))
        font = cv2.FONT_HERSHEY_SIMPLEX
        blankimage = cv2.putText(blankimage, datetimestringprint, (10,30), font,  0.9, (255,255,255), 1, cv2.LINE_AA)
        latitudestr = str('Latitude: '+str(round(Lat,3)))
        blankimage = cv2.putText(blankimage, latitudestr, (10,60), font,  0.8, (255,255,255), 1, cv2.LINE_AA)
        longitudestr = str('Longitude: '+str(round(Lon,3)))
        blankimage = cv2.putText(blankimage, longitudestr, (10,90), font,  0.8, (255,255,255), 1, cv2.LINE_AA)
        out.write(blankimage)
        cv2.imshow('Output',blankimage)
        cv2.waitKey(1)
out.release()
