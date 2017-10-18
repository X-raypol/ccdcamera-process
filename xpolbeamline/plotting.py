from astropy.visualization import (MinMaxInterval, ImageNormalize,
                                   HistEqStretch)
import matplotlib.pyplot as plt


def makeGraphs(evt):
    """Shows the different plots of the found events of the data

    Notes
    -----
    obj.analysis() must have been run before this can be run.
    The first window has positional graphs, and the second is the energy graphs
    NOTE: if the events seem to have negative KEV that means that there was to much pileup and thus when the median
    of the column is subtracted it is more than the noise, creating artificially negative events
    """

    #position graphs
    fig1 = plt.figure(figsize=(24, 18))
    fig1.canvas.set_window_title('Positional Plots of ' + self.file + self.extention)

    #makes the coordinate plot of where the events hit the detector
    main_plot = fig1.add_axes([.03, .03, .75, .6])
    im = main_plot.scatter(evt['X'], evt['Y'], c=evt['ENERGY'], marker='.',
                           edgecolor='none', cmap=plt.cm.magma)
    fig1.colorbar(im)
    main_plot.set_xlabel('X Coordinate of ROI')
    main_plot.set_ylabel('Y Coordinate of ROI')

    #makes a histogram equilized graph of the image post noise removal
    summedFiltered = self.noNoise.sum(axis=0)
    eq_plot = fig1.add_axes([.55, .65, .5, .33])
    norm = ImageNormalize(summedFiltered, interval=MinMaxInterval(),stretch=HistEqStretch(summedFiltered))
    eq_plot.imshow(summedFiltered, origin='lower', norm=norm, cmap=plt.cm.magma)

    #makes the histogram for the events along the y axis of the coord plot
    yHist = fig1.add_axes([.73, .03, .15, .6])
    yHist.hist(y,bins=335,orientation='horizontal',color='m')
    yHist.set_yticklabels([])

    #makes the histogram for the events along the x axis of the coord plot
    xHist = fig1.add_axes([.03, .65, .6, .2])
    xHist.hist(x,bins=335,color='m')
    xHist.set_xticklabels([])
    xHist.set_title('Positional Plots of Events', size=40)

    fig1.show()

    #energy graphs
    fig2 = plt.figure(figsize=(20,16))
    fig2.canvas.set_window_title('Energy analysis of '+self.file+self.extention)

    #makes the graph of the energy of the elements along the x axis
    kevPlot = fig2.add_axes([.05, .1, .65, .8])
    kevPlot.scatter(x,kev,c=kev,marker='.',edgecolor='none',cmap=plt.cm.gnuplot)
    kevPlot.set_ylabel('KEV')
    kevPlot.set_xlabel('X Coordinate of ROI')
    kevPlot.set_title('KEV vs. X of Events',size=20)

    #histogram of the the above except from 0 to 2 KEV
    biggest = np.array(kev)
    scale = int(round(biggest.max()/2*100))
    kevHist = fig2.add_axes([.75, .1, .2, .8])
    kevHist.set_xlabel('KEV Count')
    kevHist.set_title('KEV Count of Low Energy Events',size=20)
    kevHist.set_ylim([0,2])
    kevHist.hist(kev,bins=scale,orientation='horizontal',color='m')

    fig2.show()

    #self.ascaAnalysis(x,y,asca)
