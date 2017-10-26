'''This module provides a quick interface for file conversion and
quicklook display with sensible defaults for most things.
Use this while taking data on the beamline, when time is short
and a "good-enough" analysis with no knobs to turn is all you need.
'''
import os
import glob
from .sitkconverter import tiff2fitsimg
from .events import ExtractionChain
from .plotting import qualitycontrolplot

__all__ = ['UI']

class UI:
    '''Simple User Interface to the xpolbeamline package

    This class provides a simple user interface to the most important
    functionality of the xpolbeamline package.
    If provides a quick way to convert the tif files from the CCD
    camera and display a quicklook analysis to check that everything
    went well before taking the next exposure.

    An `UI` object is initialized with the path to the raw data, the
    statsfile, and an output path. After that, the methods work on the
    files in that directory.

    Parameters
    ----------
    inpath : string
        Path to the raw data (tif and txt files). Can be absolute or
        relative to the current working directory.
    statspath : string
        Path to the statfiles. This can include the filename for the
        statsfile. If the filename is not given, the converter will try
        to construct the filename using our normal file naming convention.
        Path can be absolute or relative to the current working directory.
    outpath : string or ``None``
        Path where output files are written. Can be absolute or relative
        to the current working directory.
        If not given or set to ``None``, this defaults to use the same
        location as `inpath`.
    '''
    inpath = '.'
    outpath = '.'
    statspath = '.'

    def __init__(self, inpath, statspath, outpath=None):
        self.inpath = inpath
        self.outpath = outpath if outpath is not None else inpath
        self.statspath = statspath
        self.img2evt = ExtractionChain()

    def convert(self, filename):
        '''Convert file from tif and txt to fits image and event list.

        Parameters
        ----------
        filename : string
            Name of file (no path)

        Returns
        -------
        out : string
            Filename of fits image
        evt : `astropy.table.Table`
            Events table
        '''
        print('Processing {}'.format(filename))
        out = tiff2fitsimg(filename, self.outpath, self.statspath,
                           overwrite=True)
        evt = self.img2evt(out)
        evt.write(out.replace('_img.fits', '_evt.fits'), overwrite=True)
        return out, evt

    def convert_display(self, filename):
        '''Convert tif files and display quicklook

        Parameters
        ----------
        filename : string
            Name of file (no path)

        Returns
        -------
        out : string
            Filename of fits image
        evt : `astropy.table.Table`
            Events table
        fig : `matplotlib.figure.Figure`
            Figure instance that can be used to save to file or to further modify
            the figure.
        '''
        out, evt = self.convert(filename)
        fig = qualitycontrolplot(self.img2evt.bkgremoved[0, :, :], evt,
                                 os.path.basename(out))

        return out, evt, fig

    def convert_display_newest(self):
        '''Convert most recenly modified tif file and display quicklook

        Returns
        -------
        out : string
            Filename of fits image
        evt : `astropy.table.Table`
            Events table
        '''
        tiflist = glob.glob(os.path.join(self.inpath, '*.tif'))
        if len(tiflist) == 0:
            raise FileNotFoundError('No tif file in {}'.format(self.inpath))
        newest = max(tiflist, key=os.path.getmtime)
        return self.convert_display(newest)
