import numpy as np
from datetime import datetime
from scipy.stats import sigmaclip
# offers more options, but is about 10 x slower
# from astropy.stats import sigma_clip
from scipy.ndimage import maximum_filter
from astropy.table import Table
from astropy.io import fits
from astropy.nddata.utils import extract_array
import astropy.units as u

from . import __version__ as version

grademap = np.array([[32, 64, 128],
                     [8, 0, 16],
                     [1, 2, 4]])

acis2asca = {0: [0],  # single event
             1: [1, 4, 5, 32, 128, 33, 36, 37, 129, 132, 133,
                 160, 161, 164, 165],   # diagonal split
             2: [64, 65, 68, 69, 2, 34, 130, 162],   # vertical split
             3: [8, 12, 136, 140],   # horizontal split left
             4: [16, 17, 48, 49],   # horizontal split right
             5: [3, 6, 9, 20, 40, 96, 144, 192, 13, 21, 35, 38, 44, 52, 53, 97,
                 100, 101, 131, 134, 137, 141, 145, 163, 166, 168, 172, 176,
                 177, 193, 196, 197],   # L-shaped splits
             6: [72, 76, 104, 108, 10, 11, 138, 139, 18, 22, 50, 54, 80, 81,
                 208, 209],  # L & quads
             # 7: []  # other
             }


def identify_events(fitsimage, sigma_clip_level=5):
    """Find the events of the given image.

    This function identifies local maxima, using some cuts.

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

    # subtracts the median of each column from everthing in the column
    bkgremoved = image - np.median(image, axis=1)[:, np.newaxis]

    # most important line of code below, sigma clips the data,
    # and we will assume that
    # the numbers that are clipped out at the top are events
    sc = sigmaclip(bkgremoved, high=sigma_clip_level, low=sigma_clip_level)

    # Use maximum filter to identify local peaks
    mf = maximum_filter(bkgremoved, size=(1, 3, 3))

    # masks first three columns which always have very low DN's
    # What to do with this?
    # filtered.mask[:, :, :3] = False

    frame, x, y = ((bkgremoved > sc.upper) & (bkgremoved == mf)).nonzero()

    evt = Table([x, y, frame], names=['X', 'Y', 'FRAME'])
    evt['5X5'] = [extract_array(bkgremoved[i, :, :], (5, 5), (j, k)) for i, j, k in zip(frame, x, y)]
    evt['3X3'] = evt['5X5'].data[:, 1:-1, 1:-1]
    evt['ENERGY'] = (evt['5X5'].data.sum(axis=(1, 2)) + np.random.rand(len(evt)) - .5) * 2.2 * 3.66
    evt['GRADE'] = ((evt['3X3'].data > 10) * grademap).sum(axis=(1, 2))
    evt['ASCA'] = ACIS2ASCAGrade(evt['GRADE'])
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
    # Set default to "other" and then loop over grades
    asca = 7 * np.ones_like(acis)
    for grade in acis2asca:
        asca[np.isin(acis, acis2asca[grade])] = grade

    return asca


def hotPixelFromTxt(evt, x, y):
    '''Flag all events with coordinates matching a given hot pixel list.

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
    hotpix : list of bool
        Flag column
    '''
    if len(x) != len(y):
        raise ValueError('x and y coordinates for hot pixels must have same number of entries.')

    hotpix = [(a, b) for a, b in zip(x, y)]
    return [(a, b) in hotpix for a, b in zip(evt['X'], evt['Y'])]


def hotPixelByOccurence(evt, n=None):
    """Mark events with coordinates that appear more than a threshhold number of times.

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
    hotpix : np.array of bool
        Flag column
    """
    if n is None:
        if evt.meta['Frames'] >= 9:
            n = evt.meta['Frames'] / 3

    xy = np.vstack([evt['X'], evt['Y']])
    # This is the line where everything happens:
    unique, unique_inv, unique_counts = np.unique(xy, return_inverse=True,
                                                  return_counts=True, axis=1)
    return unique_counts[unique_inv] > n


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
