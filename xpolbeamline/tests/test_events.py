import os
import tempfile
import numpy as np
from .common import tpath
from ..sitkconverter import tiff2fitsimg
from ..events import ExtractionChain


def test_convert_to_evt():
    '''Convert tif to fits with and without statsfile

    This does not test in details for correctnes,
    '''
    chain = ExtractionChain()
    with tempfile.TemporaryDirectory() as tmpdirname:
        tiff2fitsimg(tpath('wslit_-1.2deg.tif'), tmpdirname,
                     statfile=tpath(''))
        evt = chain(os.path.join(tmpdirname, 'wslit_-1.2deg_img.fits'))
    # reasonable number of events
    assert len(evt) > 20000
    # extract good evt
    good = evt[evt['GOOD']]
    assert len(good) > 15000
    # reasonable energies
    assert np.percentile(good['ENERGY'], 3.) > 0.
    assert np.percentile(good['ENERGY'], 97.) < 10000.
    # masking worked

    # events are 1 indexed
    for ax in 'XY':
        assert np.min(evt[ax]) >= 1
    # X is the longer axes
    assert np.max(evt['X']) > 1335
