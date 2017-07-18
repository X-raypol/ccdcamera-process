#Quick analysis of data
#to run open the command prompt, cd to the directory containing this and the fits files
#then run start ipython with matplotlib (type: ipython --matplotlib)
#finally type %run final_process.py

import numpy as np
from astropy.io import fits
import astropy
from astropy import stats
import matplotlib.pyplot as plt
import time
        

def sigmaClip(file):
    hdulist = fits.open(file)
    data = hdulist[0].data
    hdulist.close()

    t = time.time()

    print("Sigma Clipping") #maybe look into the scipy version
    filtered = astropy.stats.sigma_clip(data, sigma = 5)
    mask = filtered.mask.nonzero()

    print("Analyzing")
    #remove noise from events
    for frame, y, x in zip(mask[0], mask[1], mask[2]):
        if x < 3:
            filtered.data[frame][y][x] = 0
        else:
            filtered.data[frame][y][x] -= np.median(filtered.data[frame][y])
    
    #remove non events
    filtered.mask = ~filtered.mask
    noNoise = filtered.filled(0.)
    filtered.mask = ~filtered.mask

    KEVthresh = .15
    DNthresh = (KEVthresh*1000)/(2.2*3.66)

    #Collect events in a list
    events = []
    for frame, y, x in zip(mask[0], mask[1], mask[2]):
        if noNoise[frame][y][x] > DNthresh and 1339>x>0 and 1299>y>0:
            regionOfInterest = noNoise[frame,y-1:y+2,x-1:x+2].flatten()
            if noNoise[frame][y][x] == max(regionOfInterest) and sum(regionOfInterest) < 1250:#!!!!!!!!!:
                total = sum(regionOfInterest)*2.2*3.66/1000 #convert to KEV
                events.append([x,y,total])

    eventlist = np.array(events).reshape((int(len(events)),3))
    print('Time for sigma-clip %.5f seconds'%(time.time()-t)+'\n')

    #remove hot pixels
    coordList = {}
    for i in range(len(eventlist)):
        if (eventlist[i][0],eventlist[i][1]) in coordList:
            coordList[(eventlist[i][0],eventlist[i][1])] += 1
        else:
            coordList[(eventlist[i][0],eventlist[i][1])] = 1

    finalEventList = []
    framethresh = filtered.shape[0]/3
    if filtered.shape[0] < 9:
        framethresh = 3
    for j in range(len(eventlist)):
        if coordList[(eventlist[j][0],eventlist[j][1])] < framethresh:
            finalEventList.append(eventlist[j])

    #showHotPixels(coordList)
    return finalEventList

def keys(item):
    return item[1]

def showHotPixels(hotPixels):
    hotPixelsList = list(hotPixels.items())
    sortedlist = sorted(hotPixelsList, key=keys)
    for i in sortedlist:
        if i[1] > 2:
            print(i[0],i[1])


def makeGraphs(events,fig):
    x = []
    y = []
    Kev = []
    for event in events:
        x.append(event[0])
        y.append(event[1])
        Kev.append(event[2])
        
    posGraphs(x,y,fig)
    energyGraphs(x,Kev,fig)

    
def posGraphs(x,y,fig):

    #add file name titles to each plot
    plt.figure(fig*2+1)
    #make this one histogram equalized
    plt.subplot(221)
    plt.plot(x,y,'ro',markersize=2)
    plt.ylabel('Y position')
    plt.xlabel('X position')
    plt.title('Positions of Events')
    
    plt.subplot(222)
    plt.hist(y,bins=335,orientation='horizontal')
    #plt.xlabel('Y position')
    #plt.title('Occurences at each Y position')
    
    plt.subplot(223)
    plt.hist(x,bins=335)
    plt.xlabel('X position')
    #plt.title('Occurences at each X position')
    
    plt.show()
    
def energyGraphs(x,kev,fig):
    plt.figure(fig*2+2)
    
    plt.subplot(211)
    plt.plot(x,kev,'ro',markersize=2)
    plt.ylabel('KEV')
    plt.xlabel('X position')
    plt.title('X Vs. KEV')
    
    biggest = np.array(kev)
    scale = int(round(biggest.max()/2*100))
    plt.subplot(212)
    plt.xlabel('KEV Count')
    plt.xlim([0,2])
    plt.hist(kev,bins=scale)
    
    plt.show()
        

def main():
    file = input("Name of Fits file(s) (including .fits), seperate multiples with a space--> ")
    files = file.split()
    
    for i in range(len(files)):
        events = sigmaClip(files[i])
        makeGraphs(events,i*2)

if __name__ == "__main__":
    main()



