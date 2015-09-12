# -*- coding: utf-8 -*-
"""
Standalone light curve display.

For now we use a stream of random numbers as a stand-in.  I copy the scrolling
behavior of middle-right plot in scrollingPlots.py PyQtGraph example.
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.graphicsItems.ScatterPlotItem import SpotItem
import numpy as np

#SpotItem.

#setup display
win = pg.GraphicsWindow()
win.setWindowTitle('Light Curve')

#Make it a plot
p1 = win.addPlot()
#p1.setDownsampling(mode='peak')
#p1.setClipToView(True)


##We need to keep track of two things:
# - The time series data
# - The indices of selected points
data = np.empty(100)
bad = []

s1 = pg.ScatterPlotItem(brush=(255,0,0), pen='w',symbol='o')
s2 = pg.ScatterPlotItem(brush=(255,0,0), pen='b',symbol='o') #bad points
l1 = pg.PlotCurveItem()
p1.addItem(s1)
p1.addItem(s2)
p1.addItem(l1)


#Stream data
ptr = 0

def newdata():
    global data,ptr
    data[ptr] = np.random.normal()
    ptr += 1
    if ptr >= data.shape[0]:
        tmp = data
        data = np.empty(data.shape[0] * 2)
        data[:tmp.shape[0]] = tmp
    update()

def update():
    global data,ptr,bad
    goodmask=np.ones(len(data), np.bool)
    goodmask[bad] = 0
    badmask = np.zeros(len(data), np.bool)
    badmask[bad] = 1
    times = np.arange(len(data))
    s1.setData(times[goodmask[:ptr]],data[goodmask[:ptr]])
    s2.setData(times[badmask[:ptr]],data[badmask[:ptr]])
    l1.setData(times[goodmask[:ptr]],data[goodmask[:ptr]])


newdata()
timer = pg.QtCore.QTimer()
timer.timeout.connect(newdata)
timer.start(1000)

# Make points change color when clicked
## Make all plots clickable
def clicked(plot, points):
    global bad
    #print("clicked points", points)
    for p in points:
        if p.pos()[0] in bad:
            #p.resetPen()
            bad.remove(p.pos()[0])
        else:
            #p.setPen('b', width=2)
            bad.append(p.pos()[0])
    update()
s1.sigClicked.connect(clicked)
s2.sigClicked.connect(clicked)


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
