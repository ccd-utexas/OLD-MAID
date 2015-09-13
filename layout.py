# -*- coding: utf-8 -*-
"""
This scripts sets an initial layout for the ProEMOnline software.  It uses the
PyQtGraph dockarea system and was designed from the dockarea.py example.

Contains:
Left column: Observing Log
Center column: Plots
Right column: Images and Process Log
Menu bar
"""


import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.console
import numpy as np
import math
import astropy.io.fits as fits

from pyqtgraph.dockarea import *

#This program operates in four stages.
#Stage 0 - Program Initialized, waiting to open SPE file.
#Stage 1 - SPE file open, stars are being selected
#Stage 2 - Online data reduction and aperture photometry/plotting is done.
#Stage 3 - End of data acquisition detected. Final data written to file.  Timestamps verified.  Log saved.  Weather/time log data saved.
# -> revert back to Stage 0.
stage=0 #start at 0

#Keep track of the current frame:
#One version that we do science on
#One version for display purposes
def newframe(fitsfile):
    """For given filename, return science and display images.
    """
    img = fits.getdata(fitsfile)[0]
    displayimg = np.copy(img)
    #replace everything above 99%tile
    #don't do calulcations on this adjusted array!!!
    imgvals = displayimg.flatten()
    img99percentile = np.percentile(imgvals,99)
    displayimg[displayimg > img99percentile] = img99percentile
    #make color
    displayimg=np.array([displayimg,displayimg,displayimg]).transpose()
    return img,displayimg

#Start with some initial example file
fitsfile = 'ProEMExample.fits' #initial file
img,displayimg = newframe(fitsfile)

#Use a function to display a new image
#Autoscaling levels optional
def displayframe(displayimg,autoscale=False):
    """Display an RBG image
    Autoscale optional.
    Return nothing.
    """
    if autoscale:
        w5.setImage(displayimg,autoRange=True,levels=[np.min(displayimg),np.max(displayimg)-1])
    else:
        w5.setImage(displayimg,autoRange=False,autoLevels=False)

#Set up a list to keep track of star positions
starpos=[]

#Open File functionality
class WithMenu(QtGui.QMainWindow):
    
    def __init__(self):
        super(WithMenu, self).__init__()
        
        self.initUI()
        
    def initUI(self):      

        #Note: Exit is protected on Mac.  This may work on Windows.
        exitAction = QtGui.QAction('Exit', self)        
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.showDialog)

        openFile = QtGui.QAction('Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open new File')
        openFile.triggered.connect(self.showDialog)
        
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('File')
        fileMenu.addAction(openFile)
        fileMenu.addAction(exitAction)

        
    def showDialog(self):

        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', 
                '/home')
            #print str(fname)
        
        img = fits.getdata(str(fname))[0]
        w5.setImage(img)

    def closeEvent(self, event):
        
        reply = QtGui.QMessageBox.question(self, 'Message',
            "Really quit?", QtGui.QMessageBox.Yes | 
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

         
app = QtGui.QApplication([])
win = WithMenu()
area = DockArea()
win.setCentralWidget(area)
win.resize(1200,600)
win.setWindowTitle('ProEM Online Data Analysis Demo')

## Create docks, place them into the window one at a time.
## Note that size arguments are only a suggestion; docks will still have to
## fill the entire dock area and obey the limits of their internal widgets.
d1 = Dock("Dock1 - Observing Log", size=(500,300))
d2 = Dock("Dock2 - Process Log", size=(500,300))
d3 = Dock("Dock3 - Fourier Transform", size=(500,400))
d4 = Dock("Dock4 (tabbed) - Smoothed", size=(500,200))
d5 = Dock("Dock5 - Image", size=(500,200))
d6 = Dock("Dock6 (tabbed) - Light Curve", size=(500,200))
d7 = Dock("Dock7 (tabbed) - Comparison Counts", size=(500,200))
d8 = Dock("Dock8 (tabbed) - Seeing", size=(500,200))
area.addDock(d1, 'left')      ## place d1 at left edge of dock area (it will fill the whole space since there are no other docks yet)
area.addDock(d2, 'right')     ## place d2 at right edge of dock area
area.addDock(d3, 'left', d2)## place d3 at the left edge of d2
area.addDock(d4, 'top',d3)     ## place d4 on top d3
area.addDock(d5, 'top',d2)  ## place d5 on top d2
area.addDock(d6, 'above', d4)  ## place d6 above d4
area.addDock(d7, 'top', d3)   
area.addDock(d8, 'above', d7)

## Add widgets into each dock

## First dock holds the Observing Log
w1 = pg.LayoutWidget()
observer = QtGui.QLabel('Observer')
target = QtGui.QLabel('Target')
filt = QtGui.QLabel('Filter')
log = QtGui.QLabel('Log')

observerEdit = QtGui.QLineEdit()
targetEdit = QtGui.QLineEdit()
filtEdit = QtGui.QComboBox()
filtEdit.addItems(["BG40","u'","g'","r'","i'","z'","Other"])
logEdit = QtGui.QTextEdit()

w1.addWidget(observer, 1, 0)
w1.addWidget(observerEdit, 1, 1)

w1.addWidget(target, 2, 0)
w1.addWidget(targetEdit, 2, 1)
        
w1.addWidget(filt, 3, 0)
w1.addWidget(filtEdit, 3, 1)

w1.addWidget(log, 4, 0)
w1.addWidget(logEdit, 4, 1, 6, 1)

d1.addWidget(w1)


## Process Log
w2 = pg.LayoutWidget()
processLog = QtGui.QTextEdit()
processLog.setReadOnly(True)
#processLog.setTextBackgroundColor(QtGui.QColor("black"))
w2.addWidget(processLog, 0, 0, 6, 1)
d2.addWidget(w2)

## Fourier Transform - Just shows random updating noise for now
w3 = pg.PlotWidget(title="Fourier Transform")
curve = w3.plot(pen='y')
data = np.random.normal(size=(10,1000))
ptr = 0
def update():
    global curve, data, ptr, w3
    curve.setData(data[ptr%10])
    if ptr == 0:
        w3.enableAutoRange('xy', False)  ## stop auto-scaling after the first data set is plotted
    ptr += 1
timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(50)
d3.addWidget(w3)

## Smoothed Light Curve
w4 = pg.PlotWidget(title="Dock 4 plot")
w4.plot(np.random.normal(size=100))
d4.addWidget(w4)

## Image
w5 = pg.ImageView()
w5.ui.roiBtn.hide()
w5.ui.normBtn.hide()
displayframe(displayimg,autoscale=True)

def click(event):
    event.accept()
    pos = event.pos()
    #check if we're marking or unmarking a star
    #if pos.
    starpos.append([pos.x(),pos.y()])
    #img[pos.x(),pos.y()]=[255,255-img[pos.x(),pos.y(),1],255-img[pos.x(),pos.y(),1]]
    #w5.setImage(img,autoRange=False)
    processLog.append("Star selected at "+str( (int(pos.x()),int(pos.y())) ))
w5.getImageItem().mouseClickEvent = click
d5.addWidget(w5)

## Light Curve
w6 = pg.PlotWidget(title="Dock 6 plot")
w6.plot(np.random.normal(size=100))
d6.addWidget(w6)

## Smoothed Light Curve
w7 = pg.PlotWidget(title="Dock 7 plot")
w7.plot(np.random.normal(size=100))
d7.addWidget(w7)

## Smoothed Light Curve
w8 = pg.PlotWidget(title="Dock 8 plot")
w8.plot(np.random.normal(size=100))
d8.addWidget(w8)

win.show()



## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
