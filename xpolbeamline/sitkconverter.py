#cd desktop\data\final_testing

#table reading working something in the header?
#addheader function
#bordercheck better?
#test hot pixels

#add total distance form 0th order
#larger scale on graphs

#make shape graph better
#get hist eq graph binned and working
#stacking up the 5x5 of all (2 kev)  to see 5 by 5 shape

#python profiler(method timer)

import os.path
import sys
import time
import random
from glob import glob
from collections import OrderedDict
from warnings import warn
import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.time import Time
from astropy.visualization import (MinMaxInterval, ImageNormalize,
                                   HistEqStretch)
from PIL import Image
import matplotlib.pyplot as plt


statfileformat = [["Date", ""],
                  ["Time", ""],
                  ["Voltage", "volts"],
                  ["Current", "milliamps"],
                  ["PolAngle", "degrees (Polarization Angle)"],
                  ["Anode", ""],
                  ["MirPos", "millimeters (Mirror Position)"],
                  ["MirRot", "degrees (Mirror Rotation)"],
                  ["ApatTran", "centimeters (Aperture Translation)"],
                  ["GratTran", "centimeters (Grating Translation)"],
                  ["GratRot", "degrees (Grating Rotation)"],
                  ["CamTran", "centimeters (Camera Translation)"],
                  ["CamVert", "centimeters (Camera Vertical)"],
                  ["DetMirRo", "degrees (Detector Mirror Rotation)"],
                  ["GratVert", "centimeters (Grating Vertical)"],
                  ["Diode1Av", "(Diode 1 Average)"],
                  ["Diode2Av", "(Diode 2 Average)"]]
'''Column names and comments on them in statsfile'''

txtfileformat = {"Exposure": "seconds",
                 "Temp": "degrees C (of detector)",
                 "ROI": "pixels (Region of Interest)",
                 "ADC": "kilohertz (Analog-Digital Conversion)",
                 "ExpTime": "seconds (Experiment Time, as in run time)"}
'''Keys and comments for them used in txtfiles'''


class TxtParseError(Exception):
    '''Error class for problems when parsing txt files from LabView'''
    pass


class StatsFileError(Exception):
    '''Error class for problems with stats files'''
    pass


class InconsistentDataException(Exception):
    '''Error class when txt file and tiff or fits data are not consistent'''
    pass


def parse_txt(filename):
    '''Parse the txt files that LabView writes with every tiff file

    Parameters
    ----------
    filename : string
        Filename (incl path) to the txt file

    Returns
    -------
    txt : dict
        Dictionary with parsed values.
    '''
    with open(filename) as f:
        txtin = f.read()

    txtin = txtin.split(", ")
    txt = OrderdDict()
    for item in txtin:
        item = item.split(':', maxsplit=1)
        if len(item) != 2:
            raise TxtParseError('Item :{} in file {} does not follow format "key: value".'.format(item, filename))
        txt[item[0]] = item[1]

    for key in txt:
        try:
            txt[key] = int(txt[key])
        except ValueError:
            try:
                txt[key] = float(txt[key])
            except ValueError:
                pass
    # Now a few entries that need special treatment
    if 'ROI' in txt:
        txt['ROI'] = [int(x) for x in txt['ROI'].split()]

    if ('Date' in txt) and ('Time' in txt):
        month, day, year = [x.strip() for x in txt.pop('Date').split('/')]
        pm = 'PM' in txt['Time']
        txt['Time'] = txt['Time'].replace('PM', '')
        txt['Time'] = txt['Time'].replace('AM', '')
        h, m, s = [x.strip() for x in txt.pop('Time').split(':')]
        if pm:
            h = int(h) + 12
        isostr = '{year}-{month}-{day}T{h}:{m}:{s}'
        txt['DATE-OBS'] = isostr.format(year=year, month=month,
                                        day=day, h=h, m=m, s=s)
    return txt


def format_txt_as_header(txt):
    '''Format a txt dict to make it easy to insert into a fits header

    Any item with a key that is longer than 8 characters will be dropped.

    Parameters
    ----------
    txt : dict
        Dictionary with parsed values.

    Returns
    -------
    out : dict
        Dictionary where values are tuples consisting of a string and a comment
    '''
    out = OrderedDict()
    for key in txt:
        if len(key) <= 8:
            out[key] = (txt[key], txtfileformat.get(key, ''))
    if 'ROI' in txt:
        roi = out.pop('ROI')[0]
        out['ROIX0'] = (roi[0], 'Region of interest x_0 (pixel)')
        out['ROIX1'] = (roi[1], 'Region of interest x_1 (pixel)')
        out['ROIY0'] = (roi[2], 'Region of interest y_0 (pixel)')
        out['ROIY1'] = (roi[3], 'Region of interest y_1 (pixel)')
    out['ORIGIN'] = 'MIT'
    out['INSTRUME'] = ('Pol beamline', 'Polarimtery beamline')
    return out


def read_stats_file(filename):
    '''Read the stats file and parse its values.

    Parameters
    ----------
    filename : string or anything that `astropy.table.Table.read` can read
        Filename (incl path) of the statsfile

    Returns
    -------
    statstab : `astropy.table.Table`
        Table with the data from the stats file. Time and data columns are
        parsed into time objects in a new column ``datetime``
    '''
    statstab = Table.read(filename, format='ascii.no_header',
                          names=[l[0] for l in statfileformat])
    for row in statfileformat:
        statstab[row[0]].meta['comment'] = row[1]
    # Add strings together into ISO format, then parse into time
    isostr = np.char.add(np.char.add(statstab['Date'], 'T'),
                         statstab['Time'])
    statstab['datetime'] = Time(isostr)
    return statstab


def summarize_stats(start, exptime, stats, maxtdiff=5 * u.s):
    '''Find entries in statsfile for a specific dataset and summarize them.

    For numeric values, this just returns the mean of all rows in the
    statsfile which overlap the exposure. Current is treated special and
    this returns also the range of currents.

    Parameters
    ----------
    start : `astropy.time.Time`
        Start time of exposure
    exptime : `astropy.units.quantity.Quantity`
        Exposure duration
    stats : string
        Path and filename to a stats file or path to a directory of statsfiles.
        In the second case, this function will look for a file following the
        naming convention ``stats_MM_DD_YY.txt`` in the given directory.
        If such a file does not exist, all files with names starting with
        ``stats_`` will be searched for overlap with the exposure.
    maxtdiff :  `astropy.units.quantity.Quantity`
        Statsfiles are required to have entries within this time of the
        exposure start and finish. This number should be set to a time
        slightly longer than the time step size in the stats files.

    Returns
    -------
    out : dict
        Dictionary where values are tuples consisting of a string and a comment
    '''
    if os.path.isdir(stats):
        try:
            date = start.datetime
            guessfile = os.path.join(stats, 'stats_{:0=2d}_{:0=2d}_{:2d}.txt'.format(date.month, date.day, date.year-2000))
            return summarize_stats(start, exptime, guessfile, maxtdiff)
        except StatsFileError:
            statsfiles = glob(os.path.join(statfile, 'stats_*'))
            for statfile in stats:
                statstab = read_stats_file(statfile)
                overlap = ((statstab['datetime'] > start) &
                           (statstab['datetime'] < (start + exptime)))
                if overlap.sum() > 0:
                    break
        else:
            raise StatsFileError('No statsfile with matching times found in {}'.format(stats))
    elif os.path.isfile(stats):
        statfile = stats
        statstab = read_stats_file(statfile)
        overlap = ((statstab['datetime'] > start) &
                   (statstab['datetime'] < (start + exptime)))
        if overlap.sum() == 0:
            raise StatsFileError('No data in {} that matches the time of the measurement: {}'.format(statfile, start))
    else:
        raise FileNotFoundError('{} not found.'.format(stats))

    print("Loading stats file {}".format(statfile))
    overlapind = overlap.nonzero()[0]
    if (start + maxtdiff) < statstab['datetime'][overlapind[0]]:
        raise StatsFileError('No entry in statsfile {} within {} of the exposure start.'.format(statfile, maxtdiff))
    if (start + exptime - maxtdiff) > statstab['datetime'][overlapind[-1]]:
        raise StatsFileError('No entry in statsfile {} within {} of the exposure end.' +
                             'Maybe the exposure just finished and the statsfile is not yet updated?' +
                             '(Try again in a few seconds.)'.format(statfile, maxtdiff))

    out = {}
    for c in statstab.colnames:
        # Select only those that are numbers, not e.g. the Date string
        col = statstab[c]
        if hasattr(col, 'dtype') and np.issubdtype(col.dtype, np.number):
            out[c] = (np.mean(col[overlap], dtype=col.dtype),
                      col.meta['comment'])
    # Now some special treatment for columns where we want more information
    cur = statstab['Current'][overlap]
    out['CurRange'] = (np.max(cur) - np.min(cur), 'milliamps (current range)')

    return out


class SITKConverter(object):
    """Make an object to be used to converter and/or view data from the X-ray detector

    Parameters
    ----------
    file : str
        Name of the fits or tiff file (including the extention) that is to be analyzed/viewed or converted

    outpath : str
        Directory where outfiles will be written.
        If not given, this defaults to the current working directory.

    Notes
    -----
    This should primarily be run using the scipts made to run it, which in turn should be run from the command
    line using the anaconda ipython inicialized with matplotlib (unless otherwise specified)"""

    hdr = astropy.io.fits.Header()

    #this is the path to the stats file, should be changed if the dropbox file ever changes
    statsFilePath = "C:\\Users\\Polarimetry\\Dropbox (MIT)\\stats\\"

    def __init__(self, filename, outpath=None):
        #find where extention (such as ".txt") starts
        self.file, self.extention = filename.rsplit(".", maxsplit=1)
        self.inpath = os.path.dirname(self.file)
        self.file = os.path.basename(self.file)
        # find the directory of the file, such as:
        # "C:\Users\Polarimetry\Documents\Experiments"
        if outpath is not None:
            self.path = outpath
        else:
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
            The instance variable representing the image."""

        if self.extention == "tif":
            imgArray = []
            #opens the tif file in memory to be read
            tifpath = os.path.join(self.inpath, self.file + ".tif")
            with Image.open(tifpath) as img:
                #seek is the pillow function that reads tif files with multiple frames
                for i in range(img.n_frames):
                    img.seek(i)
                    #add each frame as a numpy array to a list
                    imgArray.append(np.array(img))

            #make the image a global numpy array
            self.img = np.array(imgArray)
            txtfile = os.path.join(self.inpath, self.file + '.txt')
            if os.path.isfile(txtfile):
                txt = parse_txt(txtfile)
                self.hdr.update(format_txt_as_header(txt))
                self.addName = ''
            else:
                self.addName = "_NoHeader"
                warn('No corresponding txt file with exposure data found at {}'.format(txtfile))

        elif self.extention == "fits":
            hdulist = fits.open(os.path.join(self.inpath, self.file + '.' + self.extention))
            self.img = hdulist[0].data
            self.hdr = hdulist[0].header
        else:
            raise ValueError("please run this with either a tif or a fits file")

        #simple check to make sure you are getting the number of frames you think you are
        if self.hdr['Frames'] + self.hdr['ThrowOut'] != self.img.shape[0]:
             raise InconsistentDataException("There were {} frames in this experiment, we discarded {}" +
                                             'but the Sitk thought there were {} total'.format(self.hdr['Frames'],
                                                                                               self.hdr['ThrowOut'],
                                                                                               self.img.shape[0]))
        if 'DATE-OBS' in self.hdr:
            print("Working with dataset from {}".format(self.hdr['DATE-OBS']))

    def loadAdditionalHeader(self):
        """Add data from stats file to the header data"""
        self.hdr.update(summarize_stats(Time(self.hdr['DATE-OBS']),
                                        self.hdr['ExpTime'] * u.s,
                                        self.statsFilePath))
        self.addName = ""

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
        self.filtered = astropy.stats.sigma_clip(self.img[self.throwOut: ], sigma=5)
        print('-Time for Sigma-Clip %.5f seconds'%(time.time()-t))

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
                events.append([int(x+self.xROI),int(y+self.yROI),total,int(acis),int(asca),int(frame)])
                e3x3.append(regionOfInterest)
                e5x5.append(self.noNoise[frame,y-2:y+3,x-2:x+3])

        #make final event list a numpy array
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


    def hotPixelFromTxtRemover(self):
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
            bpFinal.append(row.split())

        #for each event, only keep those whos coordinates are not in the bad pixel list
        eventList,badPixIndexs = [],[]
        for i in range(len(self.events)):
            if not [self.events[i][0],self.events[i][1]] in bpFinal:
                eventList.append(self.events[i])
                badPixIndexs.append(i)

        badPixIndexs.reverse()
        #also remove the bad pixel values from the ROI arrays
        for index in badPixIndexs:
            self.e3x3.pop(index)
            self.e5x5.pop(index)

        #update the event list
        self.events = np.array(eventList)
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
        if occurenceThresh == 3 and self.filtered.shape[0] >= 9:
            #set occurenceThresh to be 1/3 of the frames
            occurenceThresh = self.filtered.shape[0] / 3

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
        if occurenceThresh == 3 and self.filtered.shape[0] >= 9:
            #set occurenceThresh to be 1/3 of the frames
            occurenceThresh = self.filtered.shape[0] / 3

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
        while(os.path.isfile(os.path.join(self.path, newFileName))):
            newFileName= "hot_Pixels_for_"+self.file+" ("+str(counter)+").txt"
            counter += 1

        #make a new file
        outfile = open(newFileName,"w")
        #write the first line to give some information
        outfile.write("Hot Pixel List for "+self.file+", made on "+time.strftime("%d/%m/%Y")+" at "+time.strftime("%H:%M:%S"))
        #for every event in coordlist, check its occurence value
        for event in list(coordList.items()):
            #don't include those which are higher than occurenceThresh
            if event[1] >= occurenceThresh:
                outfile.write(event[0][0],event[0][1])

    def makeFITSFiles(self, overwrite=False):
        """Creates an events fits file and and an image fits file for the data

        Parameters
        ----------
        overwrite : bool
            Overwrite existing files or exit with error if files exist?

        Notes
        -----
        obj.loadData(), obj.loadAdditionalHeader(), and obj.analysis() must have been run for this to run.
        adds _table to the end of the table file and _img to the end of the image file
        """
        eventTable = Table([self.events[:, 0], self.events[:, 1],
                            self.events[:, 2],
                            self.events[:, 3], self.events[:, 4],
                            self.events[:, 5],
                            self.e3x3, self.e5x5],
                           names=['X', 'Y', 'ENERGY','ACIS','GRADE','FRAME','3X3','5X5'])
        #add the header data from the stats and info files
        eventTable['X'].unit = u.pix
        eventTable['Y'].unit = u.pix
        #eventTable['ENERGY'].unit = u.eV
        for key_value in self.hdr:
            eventTable.meta[key_value] = (self.hdr[k], self.hrd.comments[k])

        #write out the data to a fits file
        eventTable.write(os.path.join(self.path,
                                      self.file + "_events" + self.addName + ".fits"),
                         format='fits', overwrite=overwrite)

        #make a primary HDU with the image
        hdu = fits.PrimaryHDU(self.img[self.throwOut:], header=self.hdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(os.path.join(self.path,
                                     self.file + "_img" + self.addName + ".fits"),
                        overwrite=overwrite)

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
                asca.append(event[4])

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




def main():
    userText = input("Name of Fits or Tiff file(s) (including .fits or .tif), seperate multiples with a space--> ")
    files = userText.split()

    for file in files:
        fid = SITKConverter(file)
        fid.loadData()
        if str(file[file.rfind("."):]) == ".tif":
            fid.loadAdditionalHeader()
        fid.analysis()
        if str(file[file.rfind("."):]) == ".tif":
            fid.makeFITSFiles()
        fid.hotPixelByOccurenceRemover()
        fid.makeGraphs()



if __name__ == "__main__":
    main()
