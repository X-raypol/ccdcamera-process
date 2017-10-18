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
    start = Time('2017-09-27T16:33:00')
    out = summarize_stats(start, 1 * u.minute, tpath('stats_09_27_17.txt'))
    assert out['Anode'][0] == 1872
    assert np.isclose(out['Current'][0], 0.010858)


def test_stats_file_read_dir():
    '''Filename follows convention'''
    start = Time('2017-09-27T16:33:00')
    out = summarize_stats(start, 1 * u.minute, tpath(''))
    assert out['Anode'][0] == 1872
    assert np.isclose(out['Current'][0], 0.010858)


def test_stats_file_read_dir_2():
    '''Filename does not follow convention'''
    start = Time('2017-06-29T09:34:22')
    out = summarize_stats(start, 1 * u.minute, tpath(''))
    assert out['Anode'][0] == 1872
    assert np.isclose(out['Voltage'][0], 5001.)


def test_no_overlap():
    start = Time('2017-09-26T15:22:33')
    with pytest.raises(StatsFileError) as e:
        summarize_stats(start, 1 * u.minute, tpath('stats_09_27_17.txt'))
    assert 'matches the time' in str(e.value)


def test_no_statsfilefound():
    '''File does not exist'''
    start = Time('2016-09-26T15:22:33')
    with pytest.raises(FileNotFoundError):
        summarize_stats(start, 1 * u.minute, tpath('FileDoesNotExist.txt'))
