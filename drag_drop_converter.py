#drag drop converter
from SITKConverter import *

def load(path):
    fileStart = path.rfind("\\")
    fid = SITKConverter(path[fileStart:])
    fid.loadData()
    fid.loadAdditionalHeader()
    fid.analysis()
    fid.makeFITSFiles()

file_paths = sys.argv[1:]
for path in file_paths:
    load(path)
