#TODO:
#add global coordinates to be written out to the fits file for DS9
#add header elements to fits table to make it openable in DS9
#possibly stack up the 5x5 of all (2 kev) to see 5 by 5 shape

import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.visualization import (MinMaxInterval, ImageNormalize, HistEqStretch)
from PIL import Image
import matplotlib.pyplot as plt
import os.path
import sys
import time
import random
import math

class SITKConverter(object):
    """Make an object to be used to converter and/or view data from the X-ray detector

    Parameters
    ----------
    file : str
        Name of the fits or tiff file (including the extention) that is to be analyzed/viewed or converted

    Notes
    -----
    This should primarily be run using the scipts made to run it, which in turn should be run from the command
    line using the anaconda ipython inicialized with matplotlib (unless otherwise specified) (see online lab notebook
    for the full procedure)"""

    def __init__(self, file):
        #find where extention (such as ".txt") starts
        exNum = file.rfind(".")
        #cuts off the extention, so we are left with only the file name, such as "slit_0cm"
        self.file = file[:exNum]
        #save the extention
        self.extention = file[exNum:]
        #find the directory of the file, such as:
        #"C:\Users\Polarimetry\Documents\Experiments"
        self.path = os.getcwd()       
        #sets up throwOut in case we are reading a fits file image (so loadAdditionalHeader isnt called)
        self.throwOut = 0
        #sets up the ROI in case we are reading a fits file image (so loadAdditionalHeader isnt called)
        self.xROI,self.yROI = 0,0
        

    def loadData(self):
        """Reads the tiff or fits file into a numpy array.

        Returns
        -------
        self.img : 2d Numpy Array
            if the file is an image file return the instance variable representing the image.
        True : boolean
            if the file is a table inicialize the instance array with the events, then return true"""
        if self.extention == ".tif":
            #opens the tif file in memory to be read
            img = Image.open(self.file+".tif")
            imgArray = []
            #seek is the pillow function that reads tif files with multiple frames
            for i in range(img.n_frames):
                img.seek(i)
                #add each frame as a numpy array to a list
                imgArray.append(np.array(img))
            #for a check a little later
            self.totFrames = img.n_frames
            #close tif file in memory
            img.close()
            #make the image a global numpy array
            self.img = np.array(imgArray)
            #set up flag for tables
            self.table = False
            return self.img
        
        elif self.extention == ".fits":
            #opens the tif file in memory to be read
            hdulist = fits.open(self.file+self.extention)
            #check to see if the fits file is a table
            try:
                table = hdulist[1].data #this is the line that would error out if it is an img file
                #read in the relevant fields
                eventCoords = [table['X'],table['Y'],table['ENERGY'],table['GRADE']]
                self.e3x3,self.e5x5=table['3X3'],table['5X5']
                #make a global frames variable for hot pixel removal
                self.frames = np.array(table['FRAME']).max()
                #set up flag for tables
                self.table = True
                events = []
                #make an event list
                for i in range(len(eventCoords[0])):
                    events.append([eventCoords[0][i],eventCoords[1][i],eventCoords[2][i],eventCoords[3][i]])
                self.events = np.array(events)
                self.table = True
            except IndexError:
                #if it is an image file read it in
                self.img = np.asarray(hdulist[0].data, dtype=float)
                #set up flag for tables
                self.table = False 
                return self.img
        else:
            print("please run this with either a tif or a fits file")
            sys.exit()

    
    def loadAdditionalHeader(self):
        """Creates the array containing the header information of the experiment.

        Returns
        -------
        self.addHeader : 2d Numpy Array
            The instance variable containing several different header tags and then their respective data and comments

        Notes
        -----
        This is only to be used when converting a tiff file made by the labview code to a fits img and table.\n
        For this to run completely two additional files are needed, the first is the text file of the same name as the tiff file
        made at the same time as the data which contains necesary header information not contained in the tiff file itself.
        The second is the stats file (following the naming convention stats_MM_DD_YY.txt) containing the time at which the experiment was run.
        The text file should be located in the same directory as this script as well as the tiff image, and the stats file can be anywhere
        as long as the path below is updated to represent its location. Currently hard coded is the path to the dropbox."""
        
        #this is the path to the stats file, should be changed if the dropbox file ever changes
        #double backslash represents 1 backslash, it is just a special character so you have to have two
        self.statsFilePath = "C:\\Users\\Polarimetry\\Dropbox (MIT)\\stats\\"
        #gets today's date
        date = time.strftime("%x").split("/")
        #This is the name of the stats file, it should follow the convention "stats_MM_DD_YY.txt"
        self.statsFileName = "stats_"+date[0]+"_"+date[1]+"_"+date[2]+".txt"

        #shell of the array to be read into the header during fits file creation
        #format of ["header tag", "value", comment-"unit (full name of tag, if necesary)"]
        self.addHdr = [["Date","",""],
                       ["Time","",""],
                       ["Frames","",""],
                       ["Exposure","","seconds"],
                       ["Temp","","degrees C (of detector)"],
                       ["ROI","","pixels (Region of Interest)"],
                       ["ADC","","kilohertz (Analog-Digital Conversion)"],
                       ["ExpTime","","seconds (Experiment Time, as in run time)"],
                       ["Voltage","","volts"],
                       ["Current","","milliamps (Current Average)"],
                       ["CurRange","","milliamps (Current Range)"],
                       ["PolAngle","","degrees (Polarization Angle)"],
                       ["Anode","",""],
                       ["MirPos","","millimeters (Mirror Position)"],
                       ["MirRot","","degrees (Mirror Rotation)"],
                       ["ApatTran","","centimeters (Apature Translation)"],
                       ["GratTran","","centimeters (Grating Translation)"],
                       ["GratRot","","degrees (Grating Rotation)"],
                       ["CamTran","","centimeters (Camera Translation)"],
                       ["CamVert","","centimeters (Camera Vertical)"],
                       ["DetMirRo","","degrees (Detector Mirror Rotation)"],
                       ["GratVert","","centimeters (Grating Verticle)"],
                       ["Diode1Av","","(Diode 1 Average)"],
                       ["Diode2Av","","(Diode 2 Average)"]]

        self.addName = ""
                       
        #read in additional header information from the text file

        #check to make sure the file is present
        if os.path.isfile(self.path+"\\"+self.file+".txt"):
            #open the file with the same name as the tif, but with a .txt extention, in memory
            infile = open(self.file+".txt","r")
            #read in the contents of the file
            hdrData = infile.read()
            #close file in memory
            infile.close()
            #split the string on the comma into a list
            hdrDataArray = hdrData.split(", ")

            #loop through this new list
            for element in hdrDataArray:
                #for each element, headerEntry is a list of the form ["tag",data]
                headerEntry = element.split(":")
                #add the data to the correct place in selfaddhdr
                if "Date" in headerEntry:
                    self.addHdr[0][1] = headerEntry[1][1:]
                elif "Time" in headerEntry:
                    hours = 0
                    if "PM" in headerEntry[3] and int(headerEntry[1]) != 12:
                        hours = 12
                    self.addHdr[1][1] = str(int(headerEntry[1])+hours)+":"+headerEntry[2]+":"+headerEntry[3][:2]
                    #startTime is necesary for reading in the stats
                    if int(headerEntry[3][:2]) == 57 or int(headerEntry[3][:2]) == 58 or int(headerEntry[3][:2]) == 59:
                        startTime1 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2])-2)
                        startTime2 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2])-1)
                        startTime3 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2]))
                    elif int(headerEntry[3][:2]) < 10 and int(headerEntry[3][:2]) > 7:
                        startTime1 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2])-2)
                        startTime2 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2])-1)
                        startTime3 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2]))
                    elif int(headerEntry[3][:2]) <= 7:
                        startTime1 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2])+2)
                        startTime2 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2])+1)
                        startTime3 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2]))
                    else:
                        startTime1 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2]))
                        startTime2 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2])+1)
                        startTime3 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2])+2)
                elif "Frames" in headerEntry:
                    self.addHdr[2][1] = int(headerEntry[1])
                    self.frames = int(headerEntry[1])
                elif "Exposure" in headerEntry:
                    self.addHdr[3][1] = float(headerEntry[1])
                elif "Temp" in headerEntry:
                    self.addHdr[4][1] = int(headerEntry[1])
                elif "ROI" in headerEntry:
                    a = headerEntry[1].split()
                    self.addHdr[5][1] = "x1:"+a[0]+", x2:"+a[1]+", y1:"+a[2]+", y2:"+a[3]
                    self.xROI,self.yROI = round(float(a[0]))-1,round(float(a[2]))-1
                elif "ADC" in headerEntry:
                    self.addHdr[6][1] = int(headerEntry[1]) 
                elif "ExpTime" in headerEntry:
                    self.addHdr[7][1] = float(headerEntry[1])
                    #stats file is written to every 3 seconds, so the number of lines that have to do with the
                    #experiment is the length of the experiment divided by 3
                    numStatsLines = round(float(headerEntry[1]))//3
                elif "ThrowOut" in headerEntry:
                    self.throwOut = int(headerEntry[1])

            #check to make sure you are getting the number of frames you think you are
            if self.frames + self.throwOut != self.totFrames:
                print("there were ",self.frames," frames in this experiment, we through out",self.throwOut,
                      "but the Sitk thought there were",self.totFrames,"total")


            #read in additional header information from the text file

            if os.path.isfile(self.statsFilePath+"\\"+self.statsFileName):
                #index in memory to the stats file
                os.chdir(self.statsFilePath)
                #open stats file in memory
                infile = open(self.statsFileName,"r")
                #read it into a list of the different lines of the file
                statsdata = infile.readlines()
                #close the file in memory
                infile.close()
                #index in memory back to the working directory
                os.chdir(self.path)

                startIndex = 0
                #loop through the lines of the stats data looking for the start time
                for i in range(len(statsdata)):
                    if (startTime1 in statsdata[i] or startTime2 in statsdata[i] or startTime3 in statsdata[i]) and startIndex == 0:
                        startIndex = i
                #endIndex is the start time + the length of the experiment
                endIndex = startIndex + numStatsLines

                #if the start index doesn't change, something is probably wrong
                if startIndex == 0:
                    self.addName = "_NoStats"
                    print("Something went wrong when reading the stats file, one of the times "+startTime1+" / "+startTime2+" / "+startTime3+" / "
                          "was not found in the stats file,\nplease check that you have labled the stats file "
                          "correctly, or simply open today's stats file to check to see if the start time is present")

                #if the end index is larger than the last line in the file, something is probably wrong
                elif endIndex > len(statsdata):
                    self.addName = "_NoStats"
                    print("Something went wrong when reading the stats file, apparently the important stats information "
                          "starts on line",startIndex,"and ends on line",endIndex,"yet there are only",len(statsdata),
                          "lines in the file.\n please check that you have labled the stats file "
                          "correctly, or simply open today's stats file to see if there is an obvious problem")
                    
                else:
                    experimentData = []
                    #cut out the block of data where the experiment took place and split up each row
                    for j in statsdata[startIndex:endIndex]:
                        experimentData.append(j.split())
                    #make final numpy array
                    statsArray = np.array(experimentData)
                    #get non changing values from the center of the data block
                    statsArrayMiddle = len(statsArray)//2

                    #add the different values from the stats file to the header array
                    self.addHdr[8][1] = float(statsArray[statsArrayMiddle][2]) #voltage shouldn't change so get a value from the middle
                    self.addHdr[9][1] = statsArray[:,3].astype(float).mean() #Current average
                    self.addHdr[10][1] = statsArray[:,3].astype(float).max()-statsArray[:,3].astype(float).min() #current range
                    #same as voltage these values should not have changed over the course of the experiment
                    for i in range(11,22):
                        self.addHdr[i][1] = float(statsArray[statsArrayMiddle][i-7])
                    self.addHdr[22][1] = statsArray[:,15].astype(float).mean() #Diode 1 Average
                    self.addHdr[23][1] = statsArray[:,16].astype(float).mean() #Diode 2 Average
            else:
                self.addName = "_NoStats"
                print("The stats file '"+self.statsFileName+"' was not found in the directory '"
                      +self.statsFilePath+"'\nplease check to make sure it is there")

        #if the header file is not there add _NoHeader to the end of the file names
        else:
            self.addName = "_NoHeader"
            print("The header file was not found in the directory of the file being converted, currently the file "+
                  self.file+".tif is in the directory "+self.path+" please make sure that the file "+self.file+".txt "
                  "is in the same directory")
            
        return self.addHdr


    def analysis(self):
        """Find the events of the given image.

        Returns
        -------
        self.events : 2d Numpy Array
            The instance variable containing the x, y, energy, acis grating, asca grating, and frame
            of each event

        self.e3x3 : 3d Numpy Array
            The array containg the 3x3 grid of each event

        self.e5x5 : 3d Numpy Array
            The array containg the 5x5 grid of each event

        Notes
        -----
        obj.loadData() must have been run for this to run.
        The event finding algorithm used is to first use a sigma clip of 5 to find elements in the data
        that don't conform to the average, second remove the noise by subtracting the median of each column from
        the rest of the column. Finally check each value clipped out to see if it is the local maximum of its surrounding
        3x3 grid, if it is, sum the grid to get the total DN and then convert it to KEV and add the event to the array."""
        
        t = time.time() #time check to see how long the sigma clip takes is taking
        print("Sigma Clipping", end="  ") #maybe look into the scipy version
        #most important line of code below, sigma clips the data, and we will assume that
        #the numbers kept are noise values, and those clipped out are events
        self.filtered = astropy.stats.sigma_clip(self.img[self.throwOut:], sigma = 5)
        print('-Time for Sigma-Clip %.5f seconds'%(time.time()-t))
        #defined for hot pixel removal
        self.frames = self.filtered.shape[0]

        u = time.time() #time check to see how long the analysis is taking
        print("Analyzing/Finding Events", end="  ")
        #masks first three columns which always have very low DN's
        self.filtered.mask[:,:,:3] = False

        #subtracts the median of each column from everthing in the column            
        self.noNoise = self.filtered.data - np.median(self.filtered.data, axis=1)[:, np.newaxis]

        #returns an array of 3 arrays containing the x, y, and dn of all of the masked elements (the events)
        mask = self.filtered.mask.nonzero()

        events,e3x3,e5x5 = [],[],[]
        #loop through the masked 'events' and check to see if any are the highest point
        for frame, y, x in zip(mask[0], mask[1], mask[2]):
            #check whether the event is on the border of the image
            if 1 < x < 1338 and 1 < y < 1298:
                regionOfInterest = self.noNoise[frame,y-1:y+2,x-1:x+2]
                flatROI = regionOfInterest.flatten()
            else:
                flatROI,regionOfInterest = np.array([-100]),np.array([-100])
            #for each event check to see if it is the largest when compared to all neighboring pixels in a 3x3
            if self.noNoise[frame][y][x] == max(flatROI):                
                #if it is the largest, sum the 3x3, convert it to KEV
                total = (sum(flatROI)+random.random()-.5)*2.2*3.66 #convert to EV
                #calculate the ASCA grading of the event
                threshShape = regionOfInterest > 10 
                acis = sum((threshShape * np.array([[32,64,128],[8,0,16],[1,2,4]])).flatten())
                asca = self.getEVshape(acis)
                #add it to the event list
                events.append([int(x+self.xROI),int(y+self.yROI),total,int(asca),int(acis),int(frame)])
                e3x3.append(regionOfInterest)
                e5x5.append(self.noNoise[frame,y-2:y+3,x-2:x+3])
                
        #make final event lists numpy arrays
        self.events = np.array(events).reshape((int(len(events)),6))
        self.e3x3 = np.array(e3x3)
        self.e5x5 = np.array(e5x5)
        print('-Time for Analysis %.5f seconds'%(time.time()-u))
        return self.events, self.e3x3, self.e5x5

    def getEVshape(self,acis):
        """Returns asca grade given an acis grade.

        Parameters
        ----------
        acis : int
            The acis grade 0-256
        
        Returns
        -------
        asca grade : int
            The asca grade 0-7"""
        
        if acis in [0]: 
            return 0 #single event
        elif acis in [1,4,5,32,128,33,36,37,129,132,133,160,161,164,165]:
            return 1 #diagonal split
        elif acis in [64,65,68,69,2,34,130,162]:
            return 2 #vertical split
        elif acis in [8,12,136,140]:
            return 3 #horizontal split left
        elif acis in [16,17,48,49]:
            return 4 #horizontal split right
        elif acis in [3,6,9,20,40,96,144,192,13,21,35,38,44,52,53,97,100,101,131,134,137,141,145,163,166,168,172,176,177,193,196,197]:
            return 5 #L-shaped splits
        elif acis in [72,76,104,108,10,11,138,139,18,22,50,54,80,81,208,209]:
            return 6 #L & quads      
        else:
            return 7 #other
        

    def hotPixelFromTxtRemover(self, bpfile = "NoneGiven"):
        """Removes all events with whose coordinates are in a hot pixel file.
        
        Returns
        -------
        self.events : 2d Numpy Array
            The instance variable containing the x, y, energy, acis grating, asca grating, and frame
            of each event

        Notes
        -----
        obj.analysis() must have been run before this can be run.
        The text file should have 1 header line at the top then have all remaining lines be
        coordinates with the format x,y"""

        if  bpfile == "NoneGiven":
            #gets bad pixel file from the user
            bpfile = input("What is the name of the bad pixel file you would like mask the data with-> ")
        #opens it in memory
        infile = open(bpfile,"r")
        #reads the lines into an array
        badPixels = infile.readlines()
        #closes file in memory
        infile.close()

        #removes the first line as it might be meta info (such as date of list)
        badPixels.pop(0)
        #splits each line on the space to make it a list [x,y]
        bpFinal = []
        for row in badPixels:
            bpFinal.append([int(row.split(",")[0]),int(row.split(",")[1][:-1])])

        #for each event, only keep those whos coordinates are not in the bad pixel list
        eventList,badPixIndexs = [],[]
        NewE3x3,NewE5x5 = [],[]
        for i in range(len(self.events)):
            #print([str(self.events[i][0])+","+str(self.events[i][1])])
            #print([self.events[i][0],self.events[i][1]],[self.events[i][0],self.events[i][1]] in bpFinal)
            if not [self.events[i][0],self.events[i][1]] in bpFinal:
                eventList.append(self.events[i])
                NewE3x3.append(self.e3x3[i])
                NewE5x5.append(self.e5x5[i])

        #update the event list
        self.events = np.array(eventList)
        self.e3x3 = np.array(NewE3x3)
        self.e5x5 = np.array(NewE5x5)
        return self.events


    def hotPixelByOccurenceRemover(self, occurenceThresh = 3):
        """Removes all event with coordinates that appear more than a threshhold number of times.

        Parameters
        ----------
        occurenceThresh : int
            Defaults to 3, and then 1/3 of the number of frames if there are more than 9
            this is the lower bound number of times a pixel has to appear in order to be removed
        
        Returns
        -------
        self.events : 2d Numpy Array
            The instance variable containing the x, y, energy, acis grating, asca grating, and frame
            of each event

        Notes
        -----
        obj.analysis() must have been run before this can be run."""
        
        #if the user did not set an occurenceThresh and there are more than 9 frames
        if occurenceThresh == 3 and self.frames >= 9:
            #set occurenceThresh to be 1/3 of the frames
            occurenceThresh = self.frames / 3

        #make a dictionary of events with coordinates paired with occurence 
        coordList = {}
        for i in range(len(self.events)):
            #if the coordinates are in the dictionary, add one to count
            if (self.events[i][0],self.events[i][1]) in coordList:
                coordList[(self.events[i][0],self.events[i][1])] += 1
            #if they are not in the dictionary add them with a count of 1
            else:
                coordList[(self.events[i][0],self.events[i][1])] = 1

        #for each event, only keep those with a count less than occurenceThresh
        eventList = []
        for j in range(len(self.events)):
            if coordList[(self.events[j][0],self.events[j][1])] < occurenceThresh:
                eventList.append(self.events[j])

        #update the event list
        self.events = np.array(eventList)
        return self.events


    def hotPixelListMaker(self, occurenceThresh = 3):
        """Creates a hot pixel list file based on the number of occurences of the events in the event list

        Parameters
        ----------
        occurenceThresh : int
            Defaults to 3, and then 1/3 of the number of frames if there are more than 9
            this is the lower bound number of times a pixel has to appear in order to be removed
        
        Notes
        -----
        obj.analysis() must have been run before this can be run.
        Makes a file named hot_Pixels_for_<file name> containing the coordinates of the pixels
        labled "hot" based on their number of occurnces in a the data file"""

        #if the user did not set an occurenceThresh and there are more than 9 frames
        if occurenceThresh == 3 and self.frames >= 9:
            #set occurenceThresh to be 1/3 of the frames
            occurenceThresh = self.frames / 3

        #make a dictionary of events with coordinates paired with occurence 
        coordList = {}
        for i in range(len(self.events)):
            #if the coordinates are in the dictionary, add one to count
            if (self.events[i][0],self.events[i][1]) in coordList:
                coordList[(self.events[i][0],self.events[i][1])] += 1
            #if they are not in the dictionary add them with a count of 1
            else:
                coordList[(self.events[i][0],self.events[i][1])] = 1
                
        #make a new file name that doens't already exist
        newFileName = "hot_Pixels_for_"+self.file+".txt"
        counter = 1
        #check if it exists, if yes add (1), then (2) etc.
        while(os.path.isfile(self.path+"\\"+newFileName)):
            newFileName= "hot_Pixels_for_"+self.file+" ("+str(counter)+").txt"
            counter += 1

        #make a new file
        outfile = open(newFileName,"w")
        #write the first line to give some information
        outfile.write("Hot Pixel List for "+self.file+", made on "+time.strftime("%d/%m/%Y")+" at "+time.strftime("%H:%M:%S"))
        #for every event in coordlist, check its occurence value

        print(coordList)
        for event in list(coordList.items()):
            #don't include those which are higher than occurenceThresh
            if event[1] >= occurenceThresh:
                outfile.write("\n"+str(int(event[0][0]))+","+str(int(event[0][1])))



    def makeFITSFiles(self):
        """Creates an events fits file and and an image fits file for the data
        
        Notes
        -----
        obj.loadData(), obj.loadAdditionalHeader(), and obj.analysis() must have been run for this to run.
        adds _table to the end of the table file and _img to the end of the image file"""

        #make a new file name for the table that doens't already exist
        newFileName= self.file+"_events"+self.addName+".fits"
        counter = 1
        #check if it exists, if yes add (1), then (2) etc.
        while(os.path.isfile(self.path+"\\"+newFileName)):
            newFileName= self.file+"_events"+self.addName+" ("+str(counter)+").fits"
            counter += 1

        #make the table itself
        eventTable = Table([self.events[:,0], self.events[:, 1], self.events[:, 2], self.events[:, 3], self.events[:, 4], self.events[:, 5], self.e3x3, self.e5x5],
                           names=['X', 'Y', 'ENERGY','GRADE','ACIS','FRAME','3X3','5X5'])
        #add the header data from the stats and info files
        eventTable['X'].unit = u.pix
        eventTable['Y'].unit = u.pix
        for key_value in self.addHdr:
            eventTable.meta[key_value[0]]=(key_value[1],key_value[2])

        #to add header elements for the table do it here
        #eventTable.meta["name"]=("data","comment")
        
        #write out the data to a fits file
        eventTable.write(newFileName, format='fits')

        #make a new file name for the image that doens't already exist
        newFileName= self.file+"_img"+self.addName+".fits"
        counter = 1
        #check if it exists, if yes add (1), then (2) etc.
        while(os.path.isfile(self.path+"\\"+newFileName)):
            newFileName= self.file+"_img"+self.addName+" ("+str(counter)+").fits"
            counter += 1

        #make a primary HDU with the image
        hdu = fits.PrimaryHDU(self.img[self.throwOut:])
        hdulist = fits.HDUList([hdu])
        #add the header data from the stats and info files
        hdr = hdulist[0].header
        for key_value in self.addHdr:
            hdr[key_value[0]] = (key_value[1],key_value[2])
        #write out the data to a fits file
        hdulist.writeto(newFileName)

        #if the the header isn't complete ask the user if they want to keep the files
        if self.addName == "_NoHeader":
            choice = input("Delete the Tiff and Txt files, the current fits files have no header (y/n) -> ")
            if choice == "Y" or choice == "y":
                #permanemtly deletes the files (they don't go to the trash)
                os.remove(self.file+".tif")
                os.remove(self.file+".txt")
        elif self.addName == "_NoStats":
            choice = input("Delete the Tiff and Txt files, the current fits files have no stats (y/n) -> ")
            if choice == "Y" or choice == "y":
                #permanemtly deletes the files (they don't go to the trash)
                os.remove(self.file+".tif")
                os.remove(self.file+".txt")
        else:
            #permanemtly deletes the files (they don't go to the trash)
            os.remove(self.file+".tif")
            os.remove(self.file+".txt")          


    def makeGraphs(self):
        """Shows the different plots of the found events of the data
        
        Notes
        -----
        obj.analysis() must have been run before this can be run.
        The first window has positional graphs, and the second is the energy graphs
        NOTE: if the events seem to have negative KEV that means that there was to much pileup and thus when the median
        of the column is subtracted it is more than the noise, creating artificially negative events"""
        #make the three arrays to be graphed
        x,y,kev,asca = [],[],[],[]
        for event in self.events:
            if event[2]/1000 < 10:
                x.append(event[0])
                y.append(event[1])
                kev.append(event[2]/1000)
                asca.append(event[3])

        #position graphs
        fig1 = plt.figure(figsize=(24,18))
        fig1.canvas.set_window_title('Positional Plots of '+self.file+self.extention)

        #makes the coordinate plot of where the events hit the detector
        main_plot = fig1.add_axes([.03, .03, .75, .6])
        im = main_plot.scatter(x,y,c=kev,marker='.',edgecolor='none', cmap=plt.cm.magma)
        fig1.colorbar(im)
        main_plot.set_xlabel('X Coordinate of ROI')
        main_plot.set_ylabel('Y Coordinate of ROI')

        #makes a histogram equilized graph of the image post noise removal
        if self.table != True:
            summedFiltered = self.noNoise.sum(axis=0)
            eq_plot = fig1.add_axes([.55, .65, .5, .33])
            norm = ImageNormalize(summedFiltered, interval=MinMaxInterval(),stretch=HistEqStretch(summedFiltered))
            eq_plot.imshow(summedFiltered, origin='lower', norm=norm, cmap=plt.cm.magma)

        #makes the histogram for the events along the y axis of the coord plot
        yHist = fig1.add_axes([.73, .03, .15, .6])
        yHist.hist(y,bins=335,orientation='horizontal',color='m')
        yHist.set_yticklabels([])

        #makes the histogram for the events along the x axis of the coord plot
        xHist = fig1.add_axes([.03, .65, .6, .2])
        xHist.hist(x,bins=335,color='m')
        xHist.set_xticklabels([])
        xHist.set_title('Positional Plots of Events', size=40)

        fig1.show()
        
        #energy graphs
        fig2 = plt.figure(figsize=(20,16))
        fig2.canvas.set_window_title('Energy analysis of '+self.file+self.extention)

        #makes the graph of the energy of the elements along the x axis
        kevPlot = fig2.add_axes([.05, .1, .65, .8])
        kevPlot.scatter(x,kev,c=kev,marker='.',edgecolor='none',cmap=plt.cm.gnuplot)
        kevPlot.set_ylabel('KEV')
        kevPlot.set_xlabel('X Coordinate of ROI')
        kevPlot.set_title('KEV vs. X of Events',size=20)
        
        #histogram of the the above except from 0 to 2 KEV
        biggest = np.array(kev)
        scale = int(round(biggest.max()/2*100))
        kevHist = fig2.add_axes([.75, .1, .2, .8])
        kevHist.set_xlabel('KEV Count')
        kevHist.set_title('KEV Count of Low Energy Events',size=20)
        kevHist.set_ylim([0,2])
        kevHist.hist(kev,bins=scale,orientation='horizontal',color='m')
        
        fig2.show()

        #self.ascaAnalysis(x,y,asca)

    def ascaAnalysis(self,x,y,asca):
        """Experimental function to view ASCA gratings by position and shape"""
        fig3 = plt.figure(figsize=(24,18))
        ascaPlot = fig3.add_axes([.05, .05, .9, .9])
        graph = ascaPlot.scatter(x,y,c=asca,marker='.',edgecolor='none',cmap=plt.cm.gnuplot)
        fig3.colorbar(graph)
        fig3.show()

        count = [0,0,0,0,0,0,0]
        for grade in asca:
            count[int(grade)-1]=count[int(grade)-1]+1
        print("asca grade counts:")
        for i in range(len(count)):
            print(i+1,"-",count[i])

        #make a table for each grade
        #h1,bins = np.historgram(t1['x'],range = [0,1023],lim=100)
        #plt.plot[bins,hist[1,:]/hist.sum(asca=1)]
        
                       
