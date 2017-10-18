import time
import numpy as np
from datetime import datetime
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.io import fits


def identify_events(fitsimage, verbose=True):
    """Find the events of the given image.

    Parameters
    ----------
    fitsimage : sting
        Path to fits image file

    Returns
    -------
    events : `astropy.table.Table`
        Event table

    Notes
    -----
    The event finding algorithm used is to first use a sigma clip of 5 to find
    elements in the data that don't conform to the average, second remove the
    noise by subtracting the median of each column from the rest of the
    column. Finally check each value clipped out to see if it is the local
    maximum of its surrounding 3x3 grid, if it is, sum the grid to get the
    total DN and then convert it to KEV and add the event to the array.
    """
    with fits.open(fitsimage) as hdulist:
        image = hdulist[0].data
        hdr = hdulist[0].header
    if verbose:
        t = time.time()
        print("Sigma Clipping", end="  ")  # maybe look into the scipy version

    # most important line of code below, sigma clips the data,
    # and we will assume that
    # the numbers kept are noise values, and those clipped out are events
    filtered = sigma_clip(image, sigma=5)
    if verbose:
        print('-Time for Sigma-Clip %.5f seconds' % (time.time() - t))
        u = time.time()
        print("Analyzing/Finding Events", end="  ")

    # masks first three columns which always have very low DN's
    filtered.mask[:, :, :3] = False

    # subtracts the median of each column from everthing in the column
    noNoise = filtered.data - np.median(filtered.data, axis=1)[:, np.newaxis]

    # returns an array of 3 arrays containing the x, y, and dn of all of the masked elements (the events)
    mask = filtered.mask.nonzero()

    events, e3x3, e5x5 = [], [], []
    # loop through the masked 'events' and check to see if any are the highest point
    for frame, y, x in zip(mask[0], mask[1], mask[2]):
        # check whether the event is on the border of the image
        if 1 < x < 1338 and 1 < y < 1298:
            regionOfInterest = self.noNoise[frame,y-1:y+2,x-1:x+2]
            flatROI = regionOfInterest.flatten()
        else:
            flatROI, regionOfInterest = np.array([-100]), np.array([-100])
        # check if event is the largest when compared to all neighboring pixels in a 3x3
        if noNoise[frame][y][x] == max(flatROI):
            # if it is the largest, sum the 3x3, convert it to keV
            total = (sum(flatROI) + random.random() - .5) * 2.2 * 3.66
            # calculate the ASCA grading of the event
            threshShape = regionOfInterest > 10
            acis = sum((threshShape * np.array([[32, 64, 128],
                                                [8, 0, 16],
                                                [1, 2, 4]])).flatten())
            asca = ACIS2ASCAGrade(acis)
            # add it to the event list
            events.append([int(x + self.xROI), int(y + self.yROI), total,
                           int(acis), int(asca), int(frame)])
            e3x3.append(regionOfInterest)
            e5x5.append(self.noNoise[frame, y - 2: y + 3, x - 2: x + 3])

    # make final event list a numpy array
    events = np.array(events).reshape((int(len(events)), 6))
    if verbose:
        print('-Time for Analysis %.5f seconds' % (time.time() - u))

    evt = Table(events, names=['X', 'Y', 'ENERGY', 'ACIS', 'GRADE', 'FRAME'])
    evt['3X3'] = np.array(e3x3)
    evt['5X5'] = np.array(e5x5)
    evt['X'].unit = u.pix
    evt['Y'].unit = u.pix
    evt['ENERGY'].unit = u.keV

    for key_value in hdr:
        evt.meta[key_value] = (hdr[k], hdr.comments[k])

    hdr['CREATOR'] = ('XPOLBEAMLINE V' + version, 'Code for event detection')
    hdr['DATE'] = (datetime.now().isoformat(), 'Date of event detection')


def ACIS2ASCAGrade(acis):
    """Returns asca grade given an acis grade.

    Parameters
    ----------
    acis : int
        The acis grade 0-256

    Returns
    -------
    asca grade : int
        The asca grade 0-7
    """

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


def hotPixelFromTxt(evt, x, y, flagcol='hotpix'):
    '''Flag all events with whose coordinates are in a hot pixel file.

    Parameters
    ----------
    evt : `astropy.table.Table`
        Events table
    x, y : list or np.array
        List of X Y coordinates of hot pixels
    flagcol : string
        Name of column in the events file where hot pixels are recorded.
        (If the column exists it will be overwritten.)

    Returns
    -------
    evt : `astropy.table.Table`
        modified events table
    '''
    if len(x) != len(y):
        raise ValueError('x and y coordinates for hot pixels must have same number of entries.')

    hotpix = [(a, b) for a, b in zip(x, y)]
    evt[flagcol] = [(a, b) in hotpix for a, b in zip(evt['X'], evt['Y'])]
    return evt


def hotPixelByOccurence(evt, flagcol='hotpix', n=None):
    """Mark event with coordinates that appear more than a threshhold number of times.

    Parameters
    ----------
    evt : `astropy.table.Table`
        Events table
    flagcol : string
        Name of column in the events file where hot pixels are recorded.
        (If the column exists it will be overwritten.)
    n : int
        Lower bound number of times a pixel has to appear in order to be removed.
        If not set, this defaults to 3, and then 1/3 of the number of frames
        if there are more than 9.

    Returns
    -------
    evt : `astropy.table.Table`
        modified events table
    """
    if n is None:
        if evt.meta['Frames'] >= 9:
            n = evt.meta['Frames'] / 3

    xy = np.vstack([evt['X'], evt['Y']])
    # This is the line where everything happens:
    unique, unique_inv, unique_counts = np.unique(xy, return_inverse=True,
                                                  return_counts=True, axis=1)
    evt[flagcol] = unique_counts[unique_inv] > n
    return evt


def make_hotpixellist(evt, n):
    '''Write hot pixel list to file

    Parameters
    ----------
    evt : `astropy.table.Table`
        Events table
    n : int
        Lower bound number of times a pixel has to appear in order to be removed.
        If not set, this defaults to 3, and then 1/3 of the number of frames
        if there are more than 9.

    Returns
    -------
    t : `astropy.table.Table`
        Table of hot pixel coordinates
    '''
    xy = np.vstack([evt['X'], evt['Y']])
    # This is the line where everything happens:
    unique, unique_counts = np.unique(xy, return_counts=True, axis=1)
    hotpix = unique[unique_counts > n]
    t = Table(hotpix.T, names=['X', 'Y'])
    t.meta = evt.meta.copy()
    t.meta['EXTNAME'] = 'HOTPIX'
    return t
