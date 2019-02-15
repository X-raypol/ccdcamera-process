import string
import numpy as np
from datetime import datetime
from collections import OrderedDict
from scipy.stats import sigmaclip
from scipy.spatial import distance
# offers more options, but is about 10 x slower
# from astropy.stats import sigma_clip
from scipy.ndimage import maximum_filter
from astropy.table import Table
from astropy.io import fits
from astropy.nddata.utils import extract_array
from astropy.stats import sigma_clip
import astropy.units as u

from . import __version__ as version


__all__ = ['translate_wcs', 'median_column_remover', 'median_clip_remover',
           'identify_evt_sigmaclip',
           'add_islands5533', 'energy_from_island',
           'acis_grade', 'asca_grade',
           'hotpixelfromtxt', 'hotpixelbyoccurence', 'make_hotpixellist',
           'dist2nextevent',
           'ExtractionChain',
           'NoEventError',
           ]


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


class NoEventError(Exception):
    '''Exception to call if no events where found'''
    pass


def translate_wcs(evt, colnames=['X', 'Y']):
    '''Translate the WCS from image to event list

    Header keywords for WCS are different for images vs. event lists.
    This function translates from image to event list.
    It does not deal with all possible keywords, just those that
    are relevant in the context of this beamline.

    It drops the third dimension of the WCS and instead adds a time
    column based on frame number.

    Parameters
    ----------
    evt : `astropy.table.Table`
        Event table with one or more image-type WCSs in the header.
    colnames : list of two strings
        Column names in the events list that the imaging WCS applies to.
    '''
    # Fits starts counting at 1
    x = evt.colnames.index(colnames[0]) + 1
    y = evt.colnames.index(colnames[1]) + 1
    for a in [''] + list(string.ascii_letters):
        if ('WCSNAME' + a) in evt.meta:
            name = evt.meta.pop('WCSNAME' + a)
            evt.meta['TWCS{}{}'.format(x, a)] = name
            evt.meta['TWCS{}{}'.format(y, a)] = name
            evt.meta.pop('WCSAXES' + a)
            for old, new in zip(['CRVAL', 'CRPIX', 'CDELT', 'CUNIT', 'CTYPE'],
                                ['TCRVL', 'TCRPX', 'TCDLT', 'TCUNI', 'TCTYP']):
                # For all alternative WCS systems, use only 4 chars
                if len(a) == 1:
                    new = new[: -1]
                for xy, ind in zip('123', [x, y, None]):
                    val = evt.meta.pop(old + xy + a)
                    # Do not reassign Axis 3
                    if ind is not None:
                        evt.meta['{}{}{}'.format(new, ind, a)] = val
    # New time column calculated fresh from header values
    evt['TIME'] = (evt['FRAME'] + evt.meta['THROWOUT'][0]) * evt.meta['FRAMETIM'][0]
    # remove NAXIS keywords - not useful for event list
    naxlist = [k for k in evt.meta if k.startswith('NAXIS')]
    for k in naxlist:
        del evt.meta[k]
    # Table is in extension and not simple
    del evt.meta['SIMPLE']


def median_column_remover(image):
    '''Remove median for each column from image

    Parameters
    ----------
    image : np.array of shape (frame, x, y)
        original image

    Returns
    -------
    bkgremoved : np.array of same shape as ``image``
        Copy of image with background removed
    '''
    return image - np.median(image, axis=1)[:, np.newaxis, :]

def median_clip_remover(image,sigma_clip_level=3,biasout='none'):
    '''Pixel by pixel bias correction using sigma clipping

    Parameters
    ----------
    image : np.array of shape (frame, x, y)
        original image
    sigma_clip_level : integer
        how many sigma away a pixel needs to be from 
        the average to be considered an outlier
    biasout : string to 
        file path where bias will be written, e.g.
        '/nfs/cxc/h2/myhome/bias.fits' or if no
        bias needs to be written out, use 'none'

    Returns
    -------
    bkgremoved : np.array of same shape as ``image``
        Copy of image with background removed
    '''
    clipped=sigma_clip(image,axis=1)
    averaged=np.mean(clipped,axis=0)
    if biasout != 'none':
      hdu=fits.PrimaryHDU(averaged)
      hdul=fits.HDUList([hdu])
      hdul.writeto(biasout)
    return image - averaged

def identify_evt_sigmaclip(image, sigma_clip_level=5, peak_size=3):
    """Find the events of the given image.

    This function identifies local maxima, using some cuts.

    Parameters
    ----------
    image : np.array of shape (frame, x, y)
        background subtracted image
    sigma_clip_level : float
        Level for sigma clipping
    peak_size : int
        Event are recognized, if they are highest pixel in an island
        of size (peacksize * peak_size).

    Returns
    -------
    events : `astropy.table.Table`
        Event table
    """
    sc = sigmaclip(image, high=sigma_clip_level, low=sigma_clip_level)

    # Use maximum filter to identify local peaks
    mf = maximum_filter(image, size=(1, peak_size, peak_size))

    # masks first three columns which always have very low DN's
    # What to do with this?
    # filtered.mask[:, :, :3] = False

    frame, x, y = ((image > sc.upper) & (image == mf)).nonzero()

    # Fits convention is to start pixel counting at 1
    evt = Table([x + 1, y + 1, frame], names=['X', 'Y', 'FRAME'])
    evt['X'].unit = u.pix
    evt['Y'].unit = u.pix

    return evt


def add_islands5533(evt, image):
    '''Extract 5x5 and 3x3 event islands from image at FRAME, X, Y pos in ``evt``

    Parameters
    ----------
    evt : `astropy.table.Table`
        Event table
    image : np.array of shape (frame, x, y)
        3d background subtracted image
    '''
    evt['5X5'] = [extract_array(image[i, :, :], (5, 5), (j, k))
                  for i, j, k in zip(evt['FRAME'], evt['X'] - 1 - (evt.meta['ROIX0'][0] - 1),
                                     evt['Y'] - 1 - (evt.meta['ROIY0'][0] - 1))]
    evt['3X3'] = evt['5X5'].data[:, 1:-1, 1:-1]


def energy_from_island(evt, islandcol='3X3'):
    '''Returns event energy based on event island

    Parameters
    ----------
    evt : `astropy.table.Table`
        Event table
    islandcol : string
        Name of column that holds the data for the pixel island

    Returns
    -------
    energy : `astropy.units.quantity.Quantity`
        Event energy
    '''
    return (evt[islandcol].data.sum(axis=(1, 2)) +
            np.random.rand(len(evt)) - .5) * 2.2 * 3.66 *  u.eV


def acis_grade(evt):
    '''Returns acis grade based on 3x3 event island

    Parameters
    ----------
    evt : `astropy.table.Table`
        Event table

    Returns
    -------
    grade : np.array
        Acis grade 0-255
    '''
    return ((evt['3X3'].data > 10) * grademap).sum(axis=(1, 2))


def asca_grade(evt, acisgrade='GRADE'):
    """Returns asca grade given an acis grade.

    Parameters
    ----------
    evt : `astropy.table.Table`
        Event table
    acisgrade : string
        Column name in ``evt`` that holds the acis grade 0-256

    Returns
    -------
    asca grade : int
        The asca grade 0-7
    """
    # Set default to "other" and then loop over grades
    asca = 7 * np.ones_like(evt[acisgrade])
    for grade in acis2asca:
        asca[np.isin(evt[acisgrade], acis2asca[grade])] = grade

    return asca


def hotpixelfromtxt(evt, x, y):
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


def hotpixelbyoccurence(evt, n=None):
    """Mark events with coordinates that appear more than a n times.

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
        n = max(evt.meta['FRAMES'][0] / 3, 3)

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


def dist2nextevent(evt):
    d = np.zeros_like(evt['X'])
    xy = np.vstack([evt['X'].data, evt['Y'].data]).T
    for f in set(evt['FRAME']):
        ind = evt['FRAME'] == f
        sq = distance.squareform(distance.pdist(xy[ind, :], 'chebyshev'))
        np.fill_diagonal(sq, np.inf)
        d[ind] = np.min(sq, axis=1)
    return d


class ExtractionChain:
    '''Event extraction chain with good default settings

    Parameters
    ----------
    fitsimage : sting
        Path to fits image file

    Returns
    -------
    events : `astropy.table.Table`
        Event table
    '''
    bkg_remover = staticmethod(median_clip_remover)
    evt_identify = staticmethod(identify_evt_sigmaclip)
    add_islands = staticmethod(add_islands5533)
    process_steps = OrderedDict([('ENERGY', energy_from_island),
                                 ('GRADE', acis_grade),
                                 ('ASCA', asca_grade),
                                 ('HOTPIX', hotpixelbyoccurence),
                                 ('ONEDGE', lambda x: ~np.isfinite(x['ENERGY'])),
                                 ('d2nextevt', dist2nextevent),
                                 ('CONFUSED', lambda x: x['d2nextevt'] < 3),
                                 ('energy_good', lambda x: (x['ENERGY'] > 100) & (x['ENERGY'] < 10000)),
                                 ('GOOD', lambda x: ~x['ONEDGE'] & ~x['HOTPIX'] & ~x['CONFUSED'] & x['energy_good']),
                                ])

    @staticmethod
    def descr(obj):
        '''Return class or function name

        Parameters
        ----------
        obj : object

        Returns
        -------
        name : string
            name of class (for objects or classes) or function
        '''
        if hasattr(obj, "__name__"):
            return obj.__name__
        else:
            return obj.__class__.__name__

    def add_header(self, evt, hdr):
        '''Add all keywords in header to meta info of evt table.

        This function also adds a few keywords of its own.

        Parameters
        ----------
        evt : `astropy.table.Table`
            Event table
        hdr : `astropy.io.fits.Header`
            Fits header
        '''
        for k in hdr:
            # Skip reserved keyword that are automatically written
            if k in ['BITPIX']:
                continue
            evt.meta[k] = (hdr[k], hdr.comments[k])

        # Translate WCS from image to Evttable


        evt.meta['EXTNAME'] = 'EVENTS'
        evt.meta['CREATOR'] = ('XPOLBEAMLINE V' + version,
                               'Code for event detection')
        evt.meta['DATE'] = (datetime.now().isoformat(),
                            'Date of event detection')
        if 'HISTORY' not in evt.meta:
            evt.meta['HISTORY'] = []
        evt.meta['HISTORY'].append('img bkg removal function: {}'.format(self.descr(self.bkg_remover)))
        evt.meta['HISTORY'].append('evt identify function: {}'.format(self.descr(self.evt_identify)))
        evt.meta['HISTORY'].append('event islands extraction function: {}'.format(self.descr(self.add_islands)))
        translate_wcs(evt)

    def correct_xy_roi(self, evt):
        '''Add offset to X,Y coordinates based on ROI'''
        evt['X'] += self.hdr['ROIX0'] - 1
        evt['Y'] += self.hdr['ROIY0'] - 1

    def __call__(self, fitsimage):
        with fits.open(fitsimage) as hdulist:
            # Swap fits to numpy axes ordering
            self.image = np.swapaxes(np.copy(hdulist[0].data), 1, 2)
            self.hdr = hdulist[0].header

        self.image = self.image[self.hdr['THROWOUT']:, :, :]
        self.bkgremoved = self.bkg_remover(self.image)
        evt = self.evt_identify(self.bkgremoved)
        if len(evt) == 0:
            raise NoEventError('No Event was found in {} and thus event processing failed.'.format(fitsimage))
        self.add_header(evt, self.hdr)
        self.correct_xy_roi(evt)
        self.add_islands(evt, self.bkgremoved)
        if 'HISTORY' not in evt.meta:
            evt.meta['HISTORY'] = []

        for colname, call in self.process_steps.items():
            if isinstance(call, (list, tuple)) and (len(call) > 1):
                func = call[0]
                kwargs = call[1]
                evt.meta['HISTORY'].append('{} made by function {} ({})'.format(colname,
                                                                                self.descr(func),
                                                                                kwargs))
            else:
                func = call
                kwargs = {}
                evt.meta['HISTORY'].append('{} made by function {}'.format(colname, self.descr(func)))
            evt[colname] = func(evt, **kwargs)

        return evt
