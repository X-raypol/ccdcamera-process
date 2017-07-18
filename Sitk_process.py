#SITK Process
#to run this make sure you have a the tiff file, the additional header text file with
#the same name as the tiff file, and then the file with the stats called stats.txt 
#next make a batch file in the same directory with the following:
#<path to python executable> <path to this file> %*
#this is only necesary once as it is simply how this has been set up to run
#Then to drag the spe files you want converted onto this batch file

from PIL import Image
import numpy as np
from astropy.io import fits
import os.path
import sys


class Converter():
    
    def __init__(self, file, path):
        self.file = file
        self.path = path
        
    def loadImg(self):
        img = Image.open(self.file)
        self.img = np.array(img)
        img.close()

    def loadAdditionalHeader(self):
        #Additional Header from text file
        infile = open(self.file[:-4]+".txt","r")
        hdrData = infile.read()
        infile.close()
        hdrDataArray = hdrData.split(", ")
        for i in range(len(hdrDataArray)):
            hdrDataArray[i] = hdrDataArray[i].split(":")
            if hdrDataArray[i][0] == "ExpTime":
                exTimeInSeconds = round(float(hdrDataArray[i][1]))   
            elif hdrDataArray[i][0] == "Time":
                hours = 0
                if "PM" in hdrDataArray[i][3]:
                    hours = 12
                startTime = str(int(hdrDataArray[i][1])+hours)+":"+str(hdrDataArray[i][2])
                stimeH,stimeM = int(hdrDataArray[i][1])+hours, int(hdrDataArray[i][2])
                hdrDataArray[i] = ["Time",str(int(hdrDataArray[i][1])+hours)+":"+hdrDataArray[i][2]+":"+hdrDataArray[i][3][:2]]
            elif hdrDataArray[i][0] == "ROI":
                a = hdrDataArray[i][1].split()
                hdrDataArray[i] = ["ROI","x1:"+a[0]+", x2:"+a[1]+", y1:"+a[2]+", y2:"+a[3]]
        self.addHdr = hdrDataArray

        extimeH,extimeM = exTimeInSeconds//3600,exTimeInSeconds//60%60
        endTime = str(stimeH+extimeH)+":"+str(stimeM+extimeM)
        
        #Add stats file parameters to the header
        infile = open("stats.txt","r")
        statsdata = infile.readlines()
        statsdata.reverse()
        infile.close()
        
        startIndex,endIndex = 0,0
        for i in range(len(statsdata)):
            if endTime in statsdata[i] and endIndex == 0:
                endIndex = i
            elif startTime in statsdata[i]:
                startIndex = i

        experimentData = []
        for j in statsdata[endIndex:startIndex]:
            experimentData.append(j.split())
        experimentData.reverse()
        stats = np.array(experimentData)

        headerTags = ["Date","Time","Voltage","Current","PolAngle","Anode","MirPos",
                      "MirRota","ApatTran","GratTran","GratRote","CamTrans",
                      "CamVert","DetMirRo","GratVert","Diode1Av","Diode2Av"]
        
        self.addHdr.append([headerTags[2],stats[0][2]])        
        self.addHdr.append([headerTags[3],"Max: "+str(stats[:,3].astype(float).max())+", Min: "+str(stats[:,3].astype(float).min())])
        for i in range(4,14):
            self.addHdr.append([headerTags[i],stats[0][i]])
        self.addHdr.append([headerTags[15],stats[:,15].astype(float).mean()])
        self.addHdr.append([headerTags[16],stats[:,16].astype(float).mean()])
        #print(self.addHdr)
        
    def makeFITSFile(self):
        newFileName= self.file[:-4]+".fits"
        counter = 1
        while(os.path.isfile(self.path[:-len(self.file)]+newFileName)):
            newFileName= self.file[:-4]+" ("+str(counter)+").fits"
            counter += 1
        
        hdu = fits.PrimaryHDU(self.img)
        hdulist = fits.HDUList([hdu])
        hdr = hdulist[0].header
        for key_value in self.addHdr:
            hdr[key_value[0]] = key_value[1]
        
        hdulist.writeto(newFileName)

def load(fname, path):
    fid = Converter(fname, path)
    fid.loadImg()
    fid.loadAdditionalHeader()
    fid.makeFITSFile()

def getFileName(fname):
    splitFile = fname.split("\\")
    return splitFile[(len(splitFile)-1)]

file_paths = sys.argv[1:]
for p in file_paths:
    load(getFileName(p),p)
