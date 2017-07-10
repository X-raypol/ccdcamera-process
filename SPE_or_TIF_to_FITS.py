
#To run this code make a batch file in the same directory as it and put this command:
#<path to python executable> <path to this file> %*
#Then to drag the spe files you want converted onto this batch file

import numpy as N
from astropy.io import fits
import time
import os.path
import sys

class File(object):

    def __init__(self, fname, path):
        self.fname = fname
        self.path = path

    def loadTIFImg(self):
        from skimage import io
        im = io.imread(self.fname)
        return im.shape

    def loadSPEImg(self):
        self.file = open(self.fname, 'rb')
        self.shapeSetup()
        img = self.readFile(4100, self.xdim * self.ydim * self.frames, N.uint16)
        return img.reshape((self.frames, self.ydim, self.xdim))

    def shapeSetup(self):
        self.xdim = N.int64(self.readFile(42, 1, N.int16)[0])
        self.ydim = N.int64(self.readFile(656, 1, N.int16)[0])
        self.frames = N.int64(self.readFile(1446, 1, N.int32)[0])

    def readFile(self, pos, size, ntype):
        self.file.seek(pos)
        return N.fromfile(self.file, ntype, size)
    

    def makeFITSFile(self, numpyArray):
        hdu = fits.PrimaryHDU(numpyArray)
        newFileName= self.fname[:-4]+".fits"
        counter = 1
        while(os.path.isfile(self.path[:-len(self.fname)]+newFileName)):
            newFileName= self.fname[:-4]+" ("+str(counter)+").fits"
            counter += 1
        hdu.writeto(newFileName)

    def close(self):
        self.file.close()

    #if the file has a date this will get it and format it,
    #but the spe files from the SITK returned 0 at the date/time positions
    def getDateTime(self):
        rawdate = self.readFile(20, 9, N.int8)
        rawtime = self.readFile(172, 6, N.int8)
        monthLst = [("Jan","01"),("Feb","02"),("Mar","03"),("Apr","04"),("May","05"),("Jun","06"),
                     ("Jul","07"),("Aug","08"),("Sep","09"),("Oct","10"),("Nov","11"),("Dec","12")]
        monthDict = dict(monthLst)
        strdate = ''
        for ch in rawdate :
            strdate += chr(ch)
        for ch in rawtime:
            strdate += chr(ch)
        self.dateTime = strdate[5:9]+"-"+monthDict[strdate[2:5]]+"-"+strdate[0:2]+"T"+\
                        strdate[9:11]+":"+strdate[11:13]+":"+strdate[13:]
        return self.dateTime, strdate
    



def load(fname, path):
    fid = File(fname, path)
    print(fname[-4:])
    if fname[-4:] == ".spe":
        img = fid.loadSPEImg()
    elif fname[-4:] == ".tif":
        img = fid.loadTIFImg()
    fid.makeFITSFile(img)
    fid.close()

def getFileName(fname):
    splitFile = fname.split("\\")
    return splitFile[(len(splitFile)-1)]

file_paths = sys.argv[1:]
for p in file_paths:
    load(getFileName(p),p)


##if __name__ == "__main__":
##    file1 = open("newfile.txt","w")
##    file = input("SPE File to be converted ->")
##    load(file)
##    print("file converted")

    #img = load("lampe_dt.spe")
    #img = load("2017 June 14 09_39_38.spe")
    #img = load(sys.argv[-1])
    
