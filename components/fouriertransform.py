# -*- coding: utf-8 -*-
"""
Standalone smoothed light curve display.

This should basically be the same as the main lightcurve display, but it doen't
need to have the same ignore functionality.  Instead it should have a smoothing
slider to decide the kernel width.
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
from astroML.fourier import FT_continuous
from scipy.interpolate import interp1d

#setup display
win = pg.GraphicsWindow()
print type(win)
win.setWindowTitle('Smoothed Light Curve')

#Make it a plot
p1 = pg.PlotItem(labels={'left': 'amplitude', 'bottom': 'freq (per frame)'})
win.addItem(p1, row=None, col=None, rowspan=1, colspan=1)

##We need to keep track of two things:
# - The time series data
# - The indices of selected points
data = np.empty(100)
bad = []

#Set up plot components
l1 = pg.PlotCurveItem()


#Add components to plot object.
p1.addItem(l1)


#Stream data
ptr = 0
def newdata():
    global data,ptr
    data[ptr] = np.random.normal()+4.*np.sin(2.*np.pi*ptr/5.)+3.*np.sin(2.*np.pi*ptr/3.)+2.*np.sin(2.*np.pi*ptr/10000)
    ptr += 1
    #Double length of array as needed to hold new data.
    if ptr >= data.shape[0]:
        tmp = data
        data = np.empty(data.shape[0] * 2)
        data[:tmp.shape[0]] = tmp
    update()

#Update display.
def update():
    global data,ptr,bad
    #Identify which points to include/exclude
    goodmask=np.ones(len(data), np.bool)
    goodmask[bad] = 0
    times = np.arange(len(data))#Placeholder for real timestamps.
    #Calculate FT
    #replace bad values with mean
    if np.sum(goodmask[:ptr]) > 1:
        tofourier = interp1d(times[goodmask[:ptr]],data[goodmask[:ptr]])
        xnew = np.arange(min(times[goodmask[:ptr]]),max(times[goodmask[:ptr]]))
        ynew = tofourier(xnew)
        f,H = FT_continuous(xnew,ynew)
        H=2*np.sqrt(H.real**2 + H.imag**2.)/len(ynew)
        l1.setData(f[len(f)/2.:],H[len(f)/2.:])

newdata()
timer = pg.QtCore.QTimer()
timer.timeout.connect(newdata)
timer.start(20)


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
