import os
import tempfile

from ..sitkconverter import tiff2fitsimg

TEST_DIR = os.path.dirname(__file__)


def tpath(filename):
    return os.path.join(TEST_DIR, 'data', filename)


def test_convert():
    with tempfile.TemporaryDirectory() as tmpdirname:
        tiff2fitsimg(tpath('wslit_-1.0deg_+1o.tif'), tmpdirname)
