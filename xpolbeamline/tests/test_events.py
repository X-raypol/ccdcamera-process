import os
import tempfile
import numpy as np
import pytest
from .common import tpath
from ..sitkconverter import tiff2fitsimg
from ..events import ExtractionChain

chain = ExtractionChain()
with tempfile.TemporaryDirectory() as tmpdirname:
    tiff2fitsimg(tpath('wslit_-1.2deg.tif'), tmpdirname,
                 statfile=tpath(''))
    evtin = chain(os.path.join(tmpdirname, 'wslit_-1.2deg_img.fits'))

@pytest.fixture
def evt():
    return evtin.copy()

def test_convert_to_evt():
    '''Test convert gives reasonable event lists'''
    chain = ExtractionChain()
    with tempfile.TemporaryDirectory() as tmpdirname:
        tiff2fitsimg(tpath('wslit_-1.2deg.tif'), tmpdirname,
                     statfile=tpath(''))
        evt = chain(os.path.join(tmpdirname, 'wslit_-1.2deg_img.fits'))
    # reasonable number of events
    assert len(evt) > 20000
    assert evt['GOOD'].sum() > 12000

def test_convert_energies_reasonable(evt):
    good = evt[evt['GOOD']]
    # reasonable energies
    assert np.percentile(good['ENERGY'], 3.) > 0.
    assert np.percentile(good['ENERGY'], 97.) < 10000.
    # masking worked

def test_event_index(evt):
    # events are 1 indexed
    for ax in 'XY':
        assert np.min(evt[ax]) >= 1
    # X is the longer axes
    assert np.max(evt['X']) > 1335

def test_event_islands_center(evt):
    '''Look at well separated events and see if they are reasonable'''
    evt3 = evt[evt['d2nextevt'] > 5]
    # pick some random events
    for i in [4, 45, 2345]:
        row = evt[i]
        assert row['3X3'].data[1, 1] == np.max(row['3X3'])
        assert row['5X5'].data[2, 2] == np.max(row['5X5'])
        assert np.all(row['3X3'] == row['5X5'][1: -1, 1: -1])
