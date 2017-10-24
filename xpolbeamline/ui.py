import os
import glob
from .sitkconverter import tiff2fitsimg
from .events import ExtractionChain
from .plotting import qualitycontrolplot


class UI:
    inpath = '.'
    outpath = '.'
    statspath = '.'

    def __init__(self, inpath, statspath, outpath=None):
        self.inpath = inpath
        self.outpath = outpath if outpath is not None else inpath
        self.statspath = statspath
        self.img2evt = ExtractionChain()

    def convert(self, filename):
        print('Processing {}'.format(filename))
        out = tiff2fitsimg(filename, self.outpath, self.statspath,
                           overwrite=True)
        evt = self.img2evt(out)
        evt.write(out.replace('_img.fits', '_evt.fits'), overwrite=True)
        return out, evt

    def convert_display(self, filename):
        out, evt = self.convert(filename)
        fig = qualitycontrolplot(self.img2evt.bkgremoved[0, :, :], evt,
                                 os.path.basename(out))

        return out, evt, fig

    def convert_display_newest(self):
        newest = max(glob.iglob(os.path.join(self.inpath, '*.tif')),
                     key=os.path.getmtime)
        return self.convert_display(newest)
