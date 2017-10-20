import os
import tempfile
from astropy.io.fits import FITSDiff
from .common import tpath
from ..sitkconverter import tiff2fitsimg, addstats2img


def test_convert():
    '''Convert tif to fits with and without statsfile'''
    with tempfile.TemporaryDirectory() as tmpdirname:
        # convert file without stats and add those in second step
        tiff2fitsimg(tpath('wslit_-1.2deg.tif'), tmpdirname)
        outfile = os.path.join(tmpdirname, 'wslit_-1.2deg_NoStats_img.fits')
        assert os.path.isfile(outfile)
        addstats2img(outfile, statfile=tpath(''), rename=False)
        assert os.path.isfile(outfile)

        # Use statsfile directly
        outfile2 = os.path.join(tmpdirname, 'wslit_-1.2deg_img.fits')
        assert not os.path.exists(outfile2)
        tiff2fitsimg(tpath('wslit_-1.2deg.tif'), tmpdirname, statfile=tpath(''))
        assert os.path.isfile(outfile2)

        # Check result is the same
        # DATE contains time of conversion and will be slightly different
        diff = FITSDiff(outfile, outfile2, ignore_keywords=['DATE'])
        assert diff.identical


def test_add_header_rename():
    '''Check that renaming (default for addstats2img) works'''
    with tempfile.TemporaryDirectory() as tmpdirname:
        # convert file without stats and add those in second step
        tiff2fitsimg(tpath('wslit_-1.2deg.tif'), tmpdirname)
        outfile = os.path.join(tmpdirname, 'wslit_-1.2deg_NoStats_img.fits')
        outfile2 = os.path.join(tmpdirname, 'wslit_-1.2deg_img.fits')
        assert os.path.isfile(outfile)
        addstats2img(outfile, statfile=tpath(''))
        assert not os.path.exists(outfile)
        assert os.path.isfile(outfile2)
