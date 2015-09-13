# -*- coding: utf-8 -*-
"""
This scripts sets an initial layout for the ProEMOnline software.  It uses the
PyQtGraph dockarea system and was designed from the dockarea.py example.

Contains:
Left column: Observing Log
Center column: Plots
Right column: Images and Process Log
Menu bar

Keaton wrote this.
"""


import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
#import pyqtgraph.console
import numpy as np
import pickle #for saving layouts
from glob import glob
#import math
import astropy.io.fits as fits
from astroML.fourier import FT_continuous
from scipy.interpolate import interp1d

from pyqtgraph.dockarea import *

#### BEGIN PROGRAM ####



#### Define all variables that everything needs access to:

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

        #SETUP THE MENUBAR!
        #Note: Exit is protected on Mac.  This may work on Windows.
        exitAction = QtGui.QAction('Exit', self)        
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.closeEvent)

        openFile = QtGui.QAction('Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open new File')
        openFile.triggered.connect(self.showDialog)
        
        saveLayout = QtGui.QAction('Save layout', self)
        saveLayout.setStatusTip('Save the current dock layout')
        saveLayout.triggered.connect(self.saveLayout)
        
        loadLayout = QtGui.QAction('Load layout', self)
        loadLayout.setStatusTip('Load a saved dock layout')
        loadLayout.triggered.connect(self.loadLayout)
        
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('File')
        fileMenu.addAction(openFile)
        fileMenu.addAction(saveLayout)
        fileMenu.addAction(loadLayout)
        #fileMenu.addAction(exitAction)

    layoutsDir = './layouts/'
    layoutsExt = '.p'
    def saveLayout(self):
        layoutName, ok = QtGui.QInputDialog.getText(self, 'Save layout', 
            'Enter name for this layout:')
        if ok:
            pickle.dump( area.saveState(), open( self.layoutsDir+layoutName+self.layoutsExt, "wb" ) )
    
    def loadLayout(self):
        layouts = glob(self.layoutsDir+'*'+self.layoutsExt)
        if len(layouts) == 0:
            _ = QtGui.QMessageBox.warning(self,'Load layout','No saved layouts found.')
        else:
            layouts = [layout[len(self.layoutsDir):-1*len(self.layoutsExt)] for layout in layouts]
            layout, ok = QtGui.QInputDialog().getItem(self,'Load layout','Select layout: ',layouts)      
            if ok:
                state = pickle.load(open(self.layoutsDir+layout+self.layoutsExt, "rb" ) )
                area.restoreState(state)
    
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
area.addDock(d1, 'right')      ## place d1 at right edge of dock area (it will fill the whole space since there are no other docks yet)
area.addDock(d4, 'left',d1)     ## place d4 to the left of d1
area.addDock(d6, 'above', d4)  ## place d6 above d4
area.addDock(d5, 'bottom',d1)  ## place d5 below d1
area.addDock(d7, 'left', d5)   ## place d7 to the left of d5
area.addDock(d8, 'above', d7)  # place d8 above d7
area.addDock(d2, 'bottom', d5)     ## place d2 below d5
area.addDock(d3, 'left', d2) ## place d3 to the left of d2



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
w3 = pg.PlotWidget(title="Fourier Transform",labels={'left': 'amplitude', 'bottom': 'freq (per frame)'})
ft = w3.plot(pen='y')
d3.addWidget(w3)




## Smoothed Light Curve
#define smoothing function
winsize = 5.#win size in frames
def smooth(flux, window_size):
    #make sure flux is longer than window_size
    if len(flux) < window_size:
        return flux
    else:
        window = np.ones(int(window_size))/float(window_size)
        return np.convolve(flux, window, 'same')
w4 = pg.PlotWidget(title="Dock 4 plot",labels={'left': 'smoothed flux', 'bottom': 'frame #'})
ss1 = pg.ScatterPlotItem(brush=(255,0,0), pen='w',symbol='o')
sl1 = pg.PlotCurveItem()
w4.addItem(ss1)
w4.addItem(sl1)
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
w6 = pg.PlotWidget(title="Dock 6 plot",labels={'left': 'rel. flux', 'bottom': 'frame #'})

##We need to keep track of two things:
# - The time series data
# - The indices of selected points
data = np.empty(100)
bad = []

#Set up plot components
s1 = pg.ScatterPlotItem(brush=(255,0,0), pen='w',symbol='o')
s2 = pg.ScatterPlotItem(brush=(255,0,0), pen='b',symbol='o') #bad points
l1 = pg.PlotCurveItem()

#Add components to plot object.
w6.addItem(s1)
w6.addItem(s2)
w6.addItem(l1)







#Stream data
ptr = 0
def newdata():
    global data,ptr
    data[ptr] = np.random.normal()+1.*np.sin(2.*np.pi*ptr/20)+0.4*np.sin(2.*np.pi*ptr/30)+0.5*np.sin(2.*np.pi*ptr/1000)
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
    badmask = np.zeros(len(data), np.bool)
    badmask[bad] = 1
    times = np.arange(len(data))#Placeholder for real timestamps.
    s1.setData(times[goodmask[:ptr]],data[goodmask[:ptr]])
    s2.setData(times[badmask[:ptr]],data[badmask[:ptr]])
    l1.setData(times[goodmask[:ptr]],data[goodmask[:ptr]])
    ss1.setData(times[goodmask[:ptr]],smooth(data[goodmask[:ptr]],winsize))
    sl1.setData(times[goodmask[:ptr]],smooth(data[goodmask[:ptr]],winsize))
    xnew = np.arange(min(times[goodmask[:ptr]]),max(times[goodmask[:ptr]]))
    if len(xnew) > 1 and len(xnew) % 2 == 0:
        tofourier = interp1d(times[goodmask[:ptr]],data[goodmask[:ptr]])
        xnew = np.arange(min(times[goodmask[:ptr]]),max(times[goodmask[:ptr]]))
        ynew = tofourier(xnew)
        f,H = FT_continuous(xnew,ynew)
        H=2*np.sqrt(H.real**2 + H.imag**2.)/len(ynew)
        ft.setData(f[len(f)/2.:],H[len(f)/2.:])



newdata()
timer = pg.QtCore.QTimer()
timer.timeout.connect(newdata)
timer.start(100)


# Make points change color when clicked
## Make all plots clickable
def clicked(plot, points):
    global bad
    #print("clicked points", points)
    for p in points:
        if p.pos()[0] in bad:
            bad.remove(p.pos()[0])
        else:
            bad.append(p.pos()[0])
    update()
    
s1.sigClicked.connect(clicked)
s2.sigClicked.connect(clicked)

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
