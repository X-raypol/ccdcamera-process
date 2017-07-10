#quick analysis of data
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
    
    return eventlist



def makeGraphs(events):
    x = []
    y = []
    Kev = []
    for event in events:
        x.append(event[0])
        y.append(event[1])
        Kev.append(event[2])
    plot(x,y,'Positions of Events','X position','Y position',False,1)
    plot(x,Kev,'X Vs. KEV','X position','KEV',False,2)
    plot(x,False,'Occurences at each X position','X position',False,True,3)
    plot(y,False,'Occurences at each Y position','Y position',False,True,4)
    plot(Kev,False,'KEV Count','KEV',False,True,5)


def plot(x,y,title,x_axis,y_axis,hist,fig):
    plt.figure(fig)
    if hist:
        plt.hist(x,bins=670)
    else:
        plt.plot(x,y,'ro',markersize=1)
        plt.ylabel(y_axis)
    plt.title(title)
    plt.xlabel(x_axis)
    plt.show()

        

def main():
    events = sigmaClip("0cm.fits")
    makeGraphs(events)
    
    #events1 = sigmaClipV2("0cm.fits")
    #print(len(events))
    #print(events[0:15])
    #makeGraphs(events1,[5,6,7,8])


main()


fig = plt.figure()
ax = fig.add_axes([.2,.15,.7,.8])
ax = fig.add_subplot(3,1,2)

ax.plot(x,y, ...)
ax1.set_xlabel('hisytoargargfa')
ax2

fig2 = plt.figure()

ax22 = fig.add_subplot(121)

fig.axes
