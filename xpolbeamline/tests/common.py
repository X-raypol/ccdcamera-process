import os

TEST_DIR = os.path.dirname(__file__)


def tpath(filename):
    return os.path.join(TEST_DIR, 'data', filename)
