import os
import pytest
import numpy as np
from astropy.time import Time
import astropy.units as u

from ..sitkconverter import summarize_stats, StatsFileError

TEST_DIR = os.path.dirname(__file__)


def tpath(filename):
    return os.path.join(TEST_DIR, 'data', filename)


def test_summarize_stats():
    start = Time('2017-09-27T11:59:00')
    out = summarize_stats(start, 1 * u.minute, tpath('stats_09_27_17.txt'))
    assert out['Anode'][0] == 1872
    assert np.isclose(out['Voltage'][0], 5000.1)


def test_no_overlap():
    start = Time('2017-09-26T15:22:33')
    with pytest.raises(StatsFileError) as e:
        summarize_stats(start, 1 * u.minute, tpath('stats_09_27_17.txt'))
    assert 'matches the time' in str(e.value)
