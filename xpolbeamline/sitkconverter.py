from datetime import datetime
import os.path
from glob import glob
from collections import OrderedDict
import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.time import Time
from PIL import Image

from . import __version__ as version

__all__ = ['TxtParseError', 'StatsFileError', 'InconsistentDataException',
           'parse_txt', 'format_txt_as_header',
           'read_stats_file', 'summarize_stats',
           'tiff2fitsimg', 'addstats2img',
           ]

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

txtfileformat = {"Exposure": "commanded exptime per frame in seconds",
                 "Temp": "degrees C (of detector)",
                 "ROI": "pixels (Region of Interest)",
                 "ADC": "kilohertz (Analog-Digital Conversion)",
                 "ExpTime": "seconds (Experiment Time, as in run time)",
                 "Frames": "number of good frames in expsoure",
                 "ThrowOut": "number of frames to be discarded"}
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
    txt = OrderedDict()
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
        # 12:05:00 AM is 5 min after midnight and in ISO would be written
        # as 00:05:00.
        # If "PM" we add 12 h, so we need ot apply the same logic here
        # just so that 12:05 PM is still 5 min after noon, even after applying the PM correction.
        if int(h) == 12:
            h = 0
        if pm:
            h = int(h) + 12
        isostr = '{year}-{month}-{day}T{h}:{m}:{s}'
        txt['DATE-OBS'] = isostr.format(year=year, month=month,
                                        day=day, h=h, m=m, s=s)
    return txt


def format_txt_as_header(txt):
    '''Format a txt dict to make it easy to insert into a fits header

    Any item with a key that is longer than 8 characters will be dropped.
    At this point, a few generic header keywords are also added.

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
        except (StatsFileError, FileNotFoundError):
            statsfiles = glob(os.path.join(stats, 'stats_*'))
            for statfile in statsfiles:
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


def tiff2fitsimg(filename, outpath, statfile=None, overwrite=False):
    '''Convert tif files and associated txt to fits image

    This function takes a fits file and the associated txt file (assumed to
    be found in the same directory and having the same filename except
    for the extension) and converts that to a standard compliant
    fits image. It is an error if the txt file is not found.

    The imaging data is saved as a 3-D datacube in the primary extension.

    Parameters
    ----------
    filename : string
        Filename (incl path if not in the current working directory) of the tif file.
    outpath : str
        Directory where the fits image is written. The fits file will have the
        same filename as the input file (except for the extension).
    statfile : string or None
        Path and filename to a stats file or path to a directory of statsfiles.
        In the second case, this function will look for a file following the
        naming convention ``stats_MM_DD_YY.txt`` in the given directory.
        If such a file does not exist, all files with names starting with
        ``stats_`` will be searched for overlap with the exposure.
        If ``statfile=None``, "_NoStats" is appended to the name of the
        output file.
    overwrite : bool
        If the output file already exists, shall it be replaced?

    Returns
    -------
    outfile : string
        path and filename of the file that was just written
    '''
    infile = os.path.basename(filename)
    hdr = astropy.io.fits.Header()
    hdr['ORIGIN'] = 'MIT'
    hdr['INSTRUME'] = ('Pol beamline', 'MIT Polarimetry beamline')
    hdr['CREATOR'] = ('XPOLBEAMLINE V' + version, 'Code for format conversion')
    hdr['FILENAME'] = (infile, 'Original file name')
    hdr['DATE'] = (datetime.now().isoformat(), 'Date of file conversion')

    imgArray = []
    #opens the tif file in memory to be read
    with Image.open(filename) as img:
        #seek is the pillow function that reads tif files with multiple frames
        for i in range(img.n_frames):
            img.seek(i)
            #add each frame as a numpy array to a list
            imgArray.append(np.array(img))

    img = np.array(imgArray)
    txtfile = filename[:-4] + '.txt'
    txt = parse_txt(txtfile)
    hdr.update(format_txt_as_header(txt))

    # Now adjust timing keywords
    hdr.rename_keyword('EXPTIME', 'TOTLTIME')
    hdr['FRAMETIM'] = (hdr['TOTLTIME'] / (hdr['FRAMES'] + hdr['THROWOUT']),
                       'Time per recorded frame (seconds)')
    hdr['EXPTIME'] = (hdr['FRAMETIM'] * hdr['FRAMES'],
                      'Total exposure time (all non-throwout) frames')
    hdr['READTIME'] = (hdr['FRAMETIM'] - hdr['EXPOSURE'],
                       'read-out time incl in FRAMETIM (seconds)')

    if statfile is not None:
        hdr.update(summarize_stats(Time(hdr['DATE-OBS']),
                                        hdr['ExpTime'] * u.s,
                                        statfile))
        addwcs(hdr)
        addName = ""
    else:
        addName = '_NoStats'

    # simple check to make sure you are getting the correct number of frames
    if hdr['Frames'] + hdr['ThrowOut'] != img.shape[0]:
        raise InconsistentDataException("There were {} frames in this experiment, we discarded {}" +
                                        'but the Sitk thought there were {} total'.format(hdr['Frames'],
                                                                                          hdr['ThrowOut'],
                                                                                          img.shape[0]))
    if 'DATE-OBS' in hdr:
        print("Working with dataset from {}".format(hdr['DATE-OBS']))

    # make a primary HDU with the image
    hdu = fits.PrimaryHDU(img[hdr['ThrowOut']:], header=hdr)
    hdulist = fits.HDUList([hdu])
    outfile = os.path.join(outpath, infile[:-4] + addName + '_img.fits')
    hdulist.writeto(outfile, overwrite=overwrite)
    hdulist.close()
    return outfile


def addstats2img(filename, statfile, rename=True):
    '''Update header information in a fits image from the stats file

    Parameters
    ----------
    filename : string
        Filename (incl path if not in the current working directory)
    statfile : string or None
        Path and filename to a stats file or path to a directory of statsfiles.
        In the second case, this function will look for a file following the
        naming convention ``stats_MM_DD_YY.txt`` in the given directory.
        If such a file does not exist, all files with names starting with
        ``stats_`` will be searched for overlap with the exposure.
    rename : bool
        If true "_NoStats" is removed from the filename.
    '''
    with fits.open(filename, mode='update') as hdulist:
        hdr = hdulist[0].header
        hdr.update(summarize_stats(Time(hdr['DATE-OBS']), hdr['ExpTime'] * u.s,
                                   statfile))
        addwcs(hdr)
    if rename:
        os.rename(filename, filename.replace('_NoStats', ''))


def addwcs(hdr):

    hdr['WCSNAME'] = 'CAMCORD'
    hdr['WCSAXES'] = 3
    hdr['CRVAL1'] = -hdr['CAMTRAN']
    hdr['CRPIX1'] = 50
    hdr['CDELT1'] = 0.002
    hdr['CUNIT1'] = 'cm'
    hdr['CTYPE1'] = 'position'

    hdr['WCSNAMEA'] = 'DISP'
    hdr['WCSAXESA'] = 3
    hdr['CRVAL1A'] = -hdr['CAMTRAN'] * 0.15
    hdr['CRPIX1A'] = 50
    hdr['CDELT1A'] = 1.5 * hdr['CDELT1'] / 10  # 1.5 Ang / mm * pixelsize in cm / 10
    hdr['CUNIT1A'] = 'Angstroem'
    hdr['CTYPE1A'] = 'WAVE'

    for s in ['', 'A']:
        hdr['CRVAL2' + s] = hdr['CAMVERT']
        hdr['CRPIX2' + s] = 650
        hdr['CDELT2' + s] = 0.002
        hdr['CUNIT2' + s] = 'cm'
        hdr['CTYPE2' + s] = 'POSITION'
        hdr['CRVAL3' + s] = hdr['FRAMETIM'] * hdr['THROWOUT']
        hdr['CRPIX3' + s] = 1
        hdr['CDELT3' + s] = hdr['FRAMETIM']
        hdr['CUNIT3' + s] = 's'
        hdr['CTYPE3' + s] = 'TIME'
