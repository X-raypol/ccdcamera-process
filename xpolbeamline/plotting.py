import numpy as np
from astropy.visualization import (MinMaxInterval, ImageNormalize,
                                   HistEqStretch)
import matplotlib.pyplot as plt


def qualitycontrolplot(img, bkg, evt, title=''):
    fig = plt.figure(figsize=(24, 18))
    fig.canvas.set_window_title(title)

    aximg = fig.add_axes([.05, .73, .2, .22])
    norm = ImageNormalize(img, interval=MinMaxInterval(),
                          stretch=HistEqStretch(img))
    out = aximg.imshow(img.T, origin='lower', norm=norm, cmap=plt.cm.magma)
    # cbar of hist equalized plot does not work well
    # cbar = plt.colorbar(out, ax=aximg)
    aximg.set_title('Image')

    axbkg = fig.add_axes([.3, .73, .2, .22])
    norm = ImageNormalize(bkg.T, interval=MinMaxInterval())
    out = axbkg.imshow(bkg, origin='lower', norm=norm, cmap=plt.cm.magma)
    cbar = plt.colorbar(out, ax=axbkg)
    axbkg.set_title('Background')

    axasca = fig.add_axes([.65, .71, .3, .22])
    axasca.hist(evt['ASCA'], bins=np.arange(-0.5, 7.6, 1.),
                histtype='stepfilled')
    axasca.hist(evt['ASCA'][evt['GOOD']], bins=np.arange(-0.5, 7.6, 1.),
                histtype='stepfilled')
    axasca.set_title('ASCA grades')

    x0 = 0.08
    y0 = 0.08
    dx = 0.32
    dy = 0.45
    dx_s = 0.08
    dy_s = 0.15
    gap = dx + dx_s + 0.1

    axxy = fig.add_axes([x0, y0, dx, dy])
    axxyxhist = fig.add_axes([x0, y0 + dy, dx, dy_s], sharex=axxy)
    axxyyhist = fig.add_axes([x0 + dx, y0, dx_s, dy], sharey=axxy)
    plt.setp(axxyxhist.get_xticklabels(), visible=False)
    plt.setp(axxyyhist.get_yticklabels(), visible=False)

    for ind in [~evt['GOOD'], evt['GOOD']]:
        axxy.plot(evt['X'][ind], evt['Y'][ind], '.')

    xbinning = np.arange(0, 1340, 10)
    ybinning = np.arange(0, 1300, 10)
    n, bins, p = axxyxhist.hist(evt['X'], bins=xbinning, histtype='stepfilled')
    axxyxhist.hist(evt['X'][evt['GOOD']], bins=bins, histtype='stepfilled')
    n, bins, p = axxyyhist.hist(evt['Y'], bins=ybinning,
                                orientation='horizontal',
                                histtype='stepfilled', label='All events')
    axxyyhist.hist(evt['Y'][evt['GOOD']], bins=bins, orientation='horizontal',
                   histtype='stepfilled', label='good events')
    axxyyhist.legend(loc=(.1, 1.1))
    axxy.set_xlabel('X position [pixel]')
    axxy.set_ylabel('Y position [pixel]')

    axevt = fig.add_axes([x0 + gap, y0, dx, dy], sharex=axxy)
    axxhist = fig.add_axes([x0 + gap, y0 + dy, dx, dy_s], sharex=axevt)
    axyhist = fig.add_axes([x0 + dx + gap, y0, dx_s, dy], sharey=axevt)
    plt.setp(axxhist.get_xticklabels(), visible=False)
    plt.setp(axyhist.get_yticklabels(), visible=False)

    for ind in [~evt['GOOD'], evt['GOOD']]:
        axevt.plot(evt['X'][ind], evt['ENERGY'][ind] / 1e3, '.')

    n, bins, p = axxhist.hist(evt['X'], bins=xbinning, histtype='stepfilled')
    axxhist.hist(evt['X'][evt['GOOD']], bins=bins, histtype='stepfilled')
    eng = evt['ENERGY'] / 1e3
    engmin = np.nanmin(eng)
    engmax = np.nanmax(eng)
    n, bins, p = axyhist.hist(eng, bins=np.arange(engmin, engmax, .025),
                              orientation='horizontal',
                              range=[engmin, engmax], histtype='stepfilled',
                              label='All events')
    axyhist.hist(eng[evt['GOOD']], bins=bins,
                 orientation='horizontal',
                 range=[engmin, engmax], histtype='stepfilled',
                 label='good events')
    axevt.set_ylim(np.percentile(eng[evt['GOOD']], [0., 98.]))
    axevt.set_xlabel('X position [pixel]')
    axevt.set_ylabel('energy [keV]')
