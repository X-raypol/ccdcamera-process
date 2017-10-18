import argparse

from SITKConverter import SITKConverter

parser = argparse.ArgumentParser(description='''Add header info to SITK images.

For tif files, this will read both the txt file and the stats file,
for fits images, it will assume that the fits header already includes
the information form the txt file (if not, this should be added manually)
and will then add information from the stats file.

If successful, this script removes the string "_NoStats" from the filename.
''')
parser.add_argument('filename',
                    help='filename (and path) to tif or fits file.\nA txt file is expected in the same directory.')
parser.add_argument('--outpath', default='.',
                    help='Path to write output files')
parser.add_argument('--statspath', default=None,
                    help='path to directory with stats files')

args = parser.parse_args()

fid = SITKConverter(args.filename, args.outpath)
if args.statspath is not None:
    fid.statsFilePath = args.statspath
fid.loadData()
fid.loadAdditionalHeader()
fid.file = fid.file.replace('_NoStats', '')
fid.analysis()
fid.makeFITSFiles()

print("Header updated in file {} and written to {}".format(filename, args.outpath))
