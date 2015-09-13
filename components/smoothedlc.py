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

#define smoothing function
winsize = 5.
def smooth(flux, window_size):
    #make sure flux is longer than window_size
    if len(flux) < window_size:
        return flux
    else:
        window = np.ones(int(window_size))/float(window_size)
        return np.convolve(flux, window, 'same')

#setup display
win = pg.GraphicsWindow()
print type(win)
win.setWindowTitle('Smoothed Light Curve')


#Make it a plot
p1 = pg.PlotItem(labels={'left': 'rel. flux', 'bottom': 'frame #'})
win.addItem(p1, row=None, col=None, rowspan=1, colspan=1)

##We need to keep track of two things:
# - The time series data
# - The indices of selected points
data = np.empty(100)
bad = []

#Set up plot components
s1 = pg.ScatterPlotItem(brush=(255,0,0), pen='w',symbol='o')
l1 = pg.PlotCurveItem()


#Add components to plot object.
p1.addItem(s1)
p1.addItem(l1)

#Add tickslider
ts = pg.TickSliderItem(title='smooth level')
for i in range(1,11):
    ts.addTick(i)
win.nextRow()
win.addItem(ts)


#Stream data
ptr = 0
def newdata():
    global data,ptr
    data[ptr] = np.random.normal()
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
    s1.setData(times[goodmask[:ptr]],smooth(data[goodmask[:ptr]],winsize))
    l1.setData(times[goodmask[:ptr]],smooth(data[goodmask[:ptr]],winsize))

newdata()
timer = pg.QtCore.QTimer()
timer.timeout.connect(newdata)
timer.start(1000)


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
