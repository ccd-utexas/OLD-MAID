# -*- coding: utf-8 -*-
"""
Standalone light curve display.

For now we use a stream of random numbers as a stand-in.  I copy the scrolling
behavior of middle-right plot in scrollingPlots.py PyQtGraph example.
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np

win = pg.GraphicsWindow()
win.setWindowTitle('Light Curve')

p1 = win.addPlot()
p1.setDownsampling(mode='peak')
p1.setClipToView(True)
curve1 = p1.plot(pen=(200,200,200), symbolBrush=(255,0,0), symbolPen='w')

data = np.empty(100)
ptr = 0
def update():
    global data,ptr
    data[ptr] = np.random.normal()
    ptr += 1
    if ptr >= data.shape[0]:
        tmp = data
        data = np.empty(data.shape[0] * 2)
        data[:tmp.shape[0]] = tmp
    curve1.setData(data[:ptr])

timer = pg.QtCore.QTimer()
timer.timeout.connect(update)
timer.start(500)

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
