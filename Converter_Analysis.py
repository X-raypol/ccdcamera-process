#Converter/Analysis for the ccd-camera SITK process
from SITKConverter import *

def main():
    file = input("Name of Fits or Tiff file (including .fits or .tif), seperate multiples with a space--> ")
    fid = SITKConverter(file)
    fid.loadData()
    if fid.table == False:
        if str(file[file.rfind("."):]) == ".tif":
            fid.loadAdditionalHeader()
        fid.analysis()
        if str(file[file.rfind("."):]) == ".tif":
            fid.makeFITSFiles()
    #fid.hotPixelListMaker()
    #fid.hotPixelFromTxtRemover(<badPixelFile.txt>)
    fid.hotPixelByOccurenceRemover()
    fid.makeGraphs()

if __name__ == "__main__":
    main()
