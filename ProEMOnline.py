# -*- coding: utf-8 -*-
"""
This scripts sets an initial layout for the ProEMOnline software.  It uses the
PyQtGraph dockarea system and was designed from the dockarea.py example.


Keaton wrote this.
"""

#Import everything you'll need
from __future__ import absolute_import, division
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pickle #for saving layouts
from functools import partial
from glob import glob
from scipy import stats
from scipy.optimize import curve_fit
from scipy.fftpack import fft,fftfreq
import pandas as pd
import os
import subprocess
import csv
import sys
import time
import datetime as dt
import dateutil.parser
from astropy.io import fits
from scipy.interpolate import interp1d
import scipy.ndimage.filters as filters
from astropy.stats import biweight_location, biweight_midvariance
from photutils import daofind
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from pyqtgraph.dockarea import *
from bs4 import BeautifulSoup
# Local modules.
import read_spe
import mainvars as mv
import spephot
from ProcessLog import ProcessLog
from ImageDisplay import ImageDisplay
from FTPlot import FTPlot



#### BEGIN PROGRAM ####


#The organization and behavoir of the program are as follows:
#This program operates in four stages.
#Stage 0 - Program Initialized, waiting to open SPE file.
#Stage 1 - SPE file open, stars are being selected
#Stage 2 - Online data reduction and aperture photometry/plotting is being done.
#Stage 3 - End of data acquisition detected. Final data written to file.  Timestamps verified.  Log saved.  Weather/time log data saved.
# -> revert back to Stage 0.
def stagechange(num):
    if num in range(4):
        processLog.log("Program stage = "+str(num),1)
        mv.stage=num
    else: processLog.log("Attempt to change stage to invalid value ("+str(num)+")",3)


#Return a string of the current time
#Commonly needed for filenames
def timestring():
    date = dt.datetime.now()
    return date.strftime('%Y%m%d_%Hh%Mm%Ss')


#Function to save a screenshot
def saveScreenshot():
    ssfilename=os.path.splitext(mv.spefile)[0]+'_'+timestring()+'.png'
    processLog.log("Writing screenshot to file "+ssfilename,2)
    p=QtGui.QPixmap.grabWidget(area)
    writeout = p.save(ssfilename, 'png')
    if not writeout: processLog.log("Saving screenshot failed!",3)


#### STAGE 0 ####
#Set up the general GUI aspects


#Set up main window with menu items
class WithMenu(QtGui.QMainWindow):
    
    def __init__(self):
        super(WithMenu, self).__init__()
        
        self.initUI()
        
    def initUI(self):      
        #SETUP THE MENUBAR!
        #Note: Exit is protected on Mac.  This works on Windows.
        exitAction = QtGui.QAction('Exit', self)        
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)
        
        #Open SPE
        openFile = QtGui.QAction('&Open SPE', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open SPE File')
        #openFile.setCheckable(True)
        openFile.triggered.connect(self.openSPE)
        
        #Run Photometry
        runPhot = QtGui.QAction('&Run Photometry', self)
        runPhot.setShortcut('Ctrl+R')
        runPhot.setStatusTip('Run Aperture Photometry on Frames')
        runPhot.triggered.connect(self.run)
        
        #Update FT
        updateFT = QtGui.QAction('&Update FT', self)
        updateFT.setShortcut('Ctrl+U')
        updateFT.setStatusTip('Update Fourier Transform with Current Light Curve')
        updateFT.triggered.connect(self.updateFTfunct)
        
        #Run Autoguider
        autoguide = QtGui.QAction('Feed to &Autoguider', self)
        autoguide.setShortcut('Ctrl+A')
        autoguide.setStatusTip('Send most recently acquired frame to Guide82')
        autoguide.triggered.connect(self.toAutoguider)
        
        #Load dark for science frames
        loadDark = QtGui.QAction('Load Darks', self)
        loadDark.setStatusTip('Open SPE Calibrations for Dark Subtracting Science Images')
        loadDark.triggered.connect(self.openDark)
        
        #Load dark for flat frames
        loadDarkForFlats = QtGui.QAction('Load Darks for Flats', self)
        loadDarkForFlats.setStatusTip('Open SPE Calibrations for Dark Subtracting Flat Images')     
        loadDarkForFlats.triggered.connect(self.openDarkForFlats)
        
        #Load flat
        loadFlat = QtGui.QAction('Load Flats', self)
        loadFlat.setStatusTip('Open SPE Calibrations for Flatfielding Science Images')
        loadFlat.triggered.connect(self.openFlat)
        
        #Restore points
        restorePoints = QtGui.QAction('Restore Points', self)
        restorePoints.setStatusTip('Return All Previously Discarded Points to the Light Curve.')
        restorePoints.triggered.connect(self.restorePts)
                
        #Save Layout
        saveLayout = QtGui.QAction('Save Layout', self)
        saveLayout.setStatusTip('Save the current dock layout')
        saveLayout.triggered.connect(self.saveLayout)
        
        #Load Layout
        loadLayout = QtGui.QAction('Load Layout', self)
        loadLayout.setStatusTip('Load a saved dock layout')
        loadLayout.triggered.connect(self.loadLayout)
        
        #changeSmoothing
        changeSmoothing = QtGui.QAction('Change Smoothing', self)
        changeSmoothing.setStatusTip('Change Light Curve Smoothing Parameters.')
        changeSmoothing.triggered.connect(self.changeSmooth)
        
        #save screenshot
        screenshot = QtGui.QAction('Save Screenshot', self)
        screenshot.setStatusTip('Save a Screenshot of the Main Window.')
        savescreenshot = partial(saveScreenshot)
        screenshot.triggered.connect(savescreenshot)
        
        
        #Menubar
        menubar = self.menuBar()
        #File Menu
        fileMenu = menubar.addMenu('File')
        fileMenu.addAction(openFile)
        fileMenu.addAction(runPhot)
        fileMenu.addAction(updateFT)
        fileMenu.addAction(autoguide)
        fileMenu.addAction(exitAction)
        #Calibrations Menu
        calibrationsMenu = menubar.addMenu('Calibrations')
        calibrationsMenu.addAction(loadDark)
        calibrationsMenu.addAction(loadDarkForFlats)
        calibrationsMenu.addAction(loadFlat)
        #Interactions menu
        interactionsMenu = menubar.addMenu('Interact')
        interactionsMenu.addAction(restorePoints)
        self.changeApertureMenu = interactionsMenu.addMenu('Select Aperture Size')
        self.changeCompStarMenu = interactionsMenu.addMenu('Select Comp Star for Division')
        interactionsMenu.addAction(changeSmoothing)
        #Layout Menu
        layoutMenu = menubar.addMenu('Layout')
        layoutMenu.addAction(saveLayout)
        layoutMenu.addAction(loadLayout)
        #Output Menu
        outputMenu = menubar.addMenu('Output')
        outputMenu.addAction(screenshot)
        

    #Functions to save and load layouts
    layoutsDir = './layouts/'
    layoutsExt = '.p'
    def saveLayout(self):
        layoutName, ok = QtGui.QInputDialog.getText(self, 'Save layout', 
            'Enter name for this layout:')
        if ok:
            #Save dict in pickle format
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
    
    #Function to open SPE files to operate on.
    def openSPE(self):
        '''
        Select a new target SPE file to work on.
        
        Open dialog box, select file, verify that it is a SPE file.
        '''
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open SPE file', 
                mv.defaultdir,filter='Data (*.spe)'))
        if fname[-4:]=='.spe':
            processLog.log("Opening file "+fname,1)
            #Set the default directory to a couple levels up from this file
            mv.rundir = os.path.dirname(fname)
            mv.defaultdir = os.path.dirname(mv.rundir)
            #set target log text as filename to start
            targetEdit.setText(os.path.basename(fname)[:-4])
            
            #This needs to trigger a major chain of events
            stage1(fname)
        else: processLog.log("Invalid file type (must be SPE).",3)
        
    #Update the FT at user's command
    def updateFTfunct(self):
        updateft(i=mv.framenum)
        
        
    def toAutoguider(self):
        if mv.spefile != '':
            processLog.log("Opening separate program to send incoming data to Guide82.",2)
            subprocess.Popen(["python",os.path.join(os.path.dirname(os.path.abspath(__file__)),'toAutoguider.py'),mv.spefile])
        else:
            processLog.log("Open SPE file first before trying to send data to Guide82.",3)
            
    def openDark(self):
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open dark file', 
                mv.defaultdir,filter='Data (*.spe *.fits)'))
        if fname:
            processLog.log("Opening dark file "+fname,1)
            mv.dark, mv.darkExp, warnings = spephot.openDark(fname)
            if warnings: processLog.log(warnings,3)
            processLog.log("Mean dark counts: "+str(np.mean(mv.dark)))
            if mv.darkExp != mv.exptime:
                processLog.log("Exp times for dark and science frames do not match!",3)
            processLog.log("Exposure time for dark: "+str(mv.darkExp)+" s")
            processframe()
            displayFrame()
    
        
    #Load Dark frames for flat calibration
    def openDarkForFlats(self):
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open SPE dark for flat calibration', 
                mv.defaultdir,filter='Data (*.spe *.fits)'))
        if fname:
            processLog.log("Opening dark file for flats "+fname,1)
            mv.darkForFlat, mv.darkForFlatExp, warnings = spephot.openDark(fname)
            if warnings: 
                processLog.log(warnings,3)
                mv.darkForFlatExp = None #This is how we flag to not write the reduced flat to fits
            processLog.log("Mean dark counts: "+str(np.mean(mv.darkForFlat)))
            processLog.log("Exposure time for dark for flats: "+str(mv.darkForFlatExp)+" s")


    #Load Flat frames
    def openFlat(self):
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open SPE flat file', 
                mv.defaultdir,filter='Data (*.spe *.fits)'))
        if fname:
            processLog.log("Opening flat file "+fname,1)
            processLog.log("The program may hang for a moment...",1)
            mv.flat, warnings = spephot.openFlat(fname,mv.darkForFlat,mv.darkForFlatExp)
            if warnings:
                processLog.log(warnings,3)
            if mv.darkForFlatExp == 0:
                processLog.log("Bias being used for flat subtraction.",1)
            processframe()
            displayFrame()
    
    
    #Restore previously "bad" points
    def restorePts(self):
        processLog.log("Deselecting "+str(len(mv.bad))+" points.")
        mv.bad=[]
        updatelcs(i=mv.framenum)
    
    #Set up aperture size menu options
    def setupApsizeMenu(self): 
        for size in mv.apsizes:
            self.changeApertureMenu.addAction(str(size)+' pixels',lambda s=size: setApSize(s))
    
    #Set up comp star selection menu options
    def addCompStarOption(self,i): 
        self.changeCompStarMenu.addAction('Comp Star #'+str(i),lambda s=i: setCompStar(s))
    
    #Change Smoothing parameters
    def changeSmooth(self):
        kernel.openKernelDialog()
    
    #Run Photometry
    def run(self):
        #Do aperture photometry on selected stars
        if mv.stage == 1:
            if len(mv.stars) == 0:
                processLog.log("No stars selected.  Select stars before running.",3)
            else:
                mv.numstars = len(mv.stars)
                #Write original coordinates and seeing to phot_coords.orig
                f = open(mv.rundir+'/phot_coords.orig', 'w')
                for j,star in enumerate(mv.stars):
                    f.write('{:.2f} {:.2f} {:.2f}\n'.format(star[0],star[1],mv.seeing[j]))
                f.close()
                mv.selectingstars=False
                stage2()
            
        
    #Confirm Quit
    def closeEvent(self, event):
        reply = QtGui.QMessageBox.question(self, 'Message',
            "Really quit?", QtGui.QMessageBox.Yes | 
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()
            
# Make the App have a window and dock area.         
app = QtGui.QApplication([])
win = WithMenu()
area = DockArea()
win.setCentralWidget(area)
win.resize(1500,800)
win.setWindowTitle('OLD MAID Software')



## Set up each of the docks (to hold the widgets)
d1 = Dock("Observing Log", size=(500,500))
d2 = Dock("Process Log", size=(500,500))
d3 = Dock("Fourier Transform", size=(500,500))
d4 = Dock("Smoothed Light Curve", size=(1000,250))
d5 = Dock("Image", size=(500,500))
d6 = Dock("Divided Light Curve", size=(1000,250))
d7 = Dock("Raw Counts", size=(500,250))
d8 = Dock("Sky Brightness", size=(1000,250))
d9 = Dock("Seeing", size=(1000,250))

#Define initial layout
area.addDock(d4, 'left')    
area.addDock(d1, 'right',d4)    
area.addDock(d6, 'above', d4)
area.addDock(d9, 'bottom', d4)  
area.addDock(d8, 'above', d9) 
area.addDock(d7, 'above', d8)  
area.addDock(d5, 'bottom',d1)  
area.addDock(d2, 'bottom', d5)    
area.addDock(d3, 'bottom', d7) 
area.moveDock(d5,'right',d3)

#Define and place widgets into the docks

## First dock holds the Observing Log
#Type of widget: Form
w1 = pg.LayoutWidget()
#Name the form elements
observer = QtGui.QLabel('Observer')
target = QtGui.QLabel('Target')
filt = QtGui.QLabel('Filter')
logtext = QtGui.QLabel('Log')
#Define the types of fields
observerEdit = QtGui.QLineEdit()
targetEdit = QtGui.QLineEdit()
filtEdit = QtGui.QComboBox()
filtEdit.addItems(["BG40","u'","g'","r'","i'","z'","Other"])
logEdit = QtGui.QTextEdit()
logEdit.setText("WARNING: None of these log fields are saved!")
#Put the fields in the form
w1.addWidget(observer, 1, 0)
w1.addWidget(observerEdit, 1, 1)
w1.addWidget(target, 2, 0)
w1.addWidget(targetEdit, 2, 1)        
w1.addWidget(filt, 3, 0)
w1.addWidget(filtEdit, 3, 1)
w1.addWidget(logtext, 4, 0)
w1.addWidget(logEdit, 4, 1, 6, 1)
#Put the widget in the dock
d1.addWidget(w1)


## Process Log
# Records activity.
w2 = pg.LayoutWidget()
processLog = ProcessLog()
w2.addWidget(processLog, 0, 0, 6, 1)
d2.addWidget(w2)

## Light Curve
# It's a plot
w6 = pg.PlotWidget(title="Divided Light Curve",labels={'left': 'rel. flux', 'bottom': 'time (s)'})
# Set up plot components
# Raw points
s1 = pg.ScatterPlotItem(brush=(255,0,0), pen='w',symbol='o')
# Bad (ignored) points #Not currently displayed since it causes scaling issues.
#s2 = pg.ScatterPlotItem(brush=(255,0,0), pen='b',symbol='o')
# Connecting lines
l1 = pg.PlotCurveItem()
#Add components to plot widget.
w6.addItem(s1)
#w6.addItem(s2)
w6.addItem(l1)
#Add widget to dock
d6.addWidget(w6)
# Make points change color when clicked
def clicked(plot, points):
    for p in points:
        if p.pos()[0]/mv.exptime in mv.bad:
            mv.bad.remove(p.pos()[0]/mv.exptime)
        else:
            mv.bad.append(p.pos()[0]/mv.exptime)
    updatelcs(i=mv.framenum)
    
s1.sigClicked.connect(clicked)
#s2.sigClicked.connect(clicked)


## Smoothed Light Curve
w4 = pg.PlotWidget(title="Smoothed Light Curve",labels={'left': 'smoothed flux', 'bottom': 'time (s)'})
ss1 = pg.ScatterPlotItem(brush=(255,0,0), pen='w',symbol='o')
sl1 = pg.PlotCurveItem()
w4.addItem(ss1)
w4.addItem(sl1)
d4.addWidget(w4)


## Raw Star/Sky Counts
w7 = pg.PlotWidget(title="Raw Star Counts",labels={'left': 'flux summed in aperture', 'bottom': 'time (s)'})
d7.addWidget(w7)
#Hold the individual plot items in this list once they are created:
rawcounts=[]

## Sky
w8 = pg.PlotWidget(title="Sky Brightness",labels={'left': 'median sky counts', 'bottom': 'time (s)'})
sky = pg.PlotCurveItem()
w8.addItem(sky)
d8.addWidget(w8)

## Seeing
w9 = pg.PlotWidget(title="Seeing",labels={'left': 'FWHM (pixels)', 'bottom': 'time (s)'})
d9.addWidget(w9)
gridlines = pg.GridItem()
w9.addItem(gridlines)
#Hold the individual plot items in this list once they are created:
seeingplots = []


## Fourier Transform
ftplot = FTPlot()
d3.addWidget(ftplot)


## Image
imageDisplay = ImageDisplay()
#w5.ui.normBtn.hide() #Causing trouble on windows
#Define function for selecting stars. (must be defined before linking the click action)
def click(event):#Linked to image click event
    if event.button() == 1 and mv.selectingstars:
        event.accept()
        pos = event.pos()
        #x and y are swapped in the GUI!
        x=pos.x()
        y=pos.y()
        #processLog.log('Clicked at ({:.2f}, {:.2f})'.format(x,y),level=0)
        #improve coordinates
        dx,dy,newseeing = improvecoords(x,y)
        #round originals so original position *within* pixel doesn't affect answer
        newcoords=[np.floor(x)+dx,np.floor(y)+dy]
        mv.stars.append(newcoords)
        mv.seeing.append(newseeing)
        #make menuoption for comp star selection
        if len(mv.stars) > 1: win.addCompStarOption(len(mv.stars)-1)
        #Mark stars in image display
        imageDisplay.displayApertures(mv.stars,4)
        #Set up plot for raw counts and seeing:
        rawcounts.append(pg.ScatterPlotItem(pen=pencolors[len(mv.stars)-1],symbol='o',size=1))
        seeingplots.append(pg.PlotCurveItem(pen=seeingcolors[len(mv.stars)-1]))
        processLog.log('Star selected at ({:.2f}, {:.2f})'.format(newcoords[0],newcoords[1]),level=1)
        
    elif event.button() == 2: 
        event.accept()#Passed on to other functionality if not accepted.
        print "RIGHT!"


imageDisplay.getImageItem().mouseClickEvent = click #Function defined below
#w5.keyPressEvent = moveCircles # Seems to be the right thing for detecting frame changes,
#But I can't connect to it without overriding other behavior.  May need to subclass this.


#Set up plot for apertures around stars
#print QtGui.QColor.colorNames() for available names.
stringcolors=['red','green','blue','magenta','orange','yellow',
              'darkred','darkgreen','darkblue','darkmagenta','darkorange','darkgoldenrod',
              'hotpink','seagreen','skyblue','salmon','brown','lightyellow']
pencolors = [pg.mkPen(QtGui.QColor(c), width=3) for c in stringcolors]
seeingcolors = [pg.mkPen(QtGui.QColor(c), width=1.5) for c in stringcolors]
#Add widget to dock

d5.addWidget(imageDisplay)



## Show the program!
win.show()
win.raise_()

#win.activateWindow()













# I think everything is set up enough to start doing stuff
# Send initial message to process log.
processLog.log("ProEMOnline initialized",2)
#log("Development version.  Do not trust.",3)
stagechange(0)
processLog.log("Open SPE file to begin analysis.",1)



#### STAGE 1 ####

# Stage 1 starts when a SPE file is loaded.
# It's the "getting everything set up" stage
# Since the SPE file is loaded by the menu action, this will be one big
# function that is called on the new image.


def stage1(fname):
    #Load SPE File
    #Announce Stage 1    
    stagechange(1)
    #Record SPE filename this once
    mv.spefile = fname
    #Read in SPE data
    mv.spe = read_spe.File(mv.spefile)
    mv.binning = 1024/mv.spe.get_frame(0)[0].shape[0]
    processLog.log(str(mv.spe.get_num_frames()) + ' frames read in.')
    mv.exptime=getexptime(mv.spe)
    processLog.log('Inferred exposure time: '+str(mv.exptime)+' s')
    if hasattr(mv.spe, 'footer_metadata'): 
        #log('SPE file has footer.')
        mv.exptime=np.round(float(BeautifulSoup(mv.spe.footer_metadata, "xml").find(name='ExposureTime').text)/1000.)
        #log('Exposute time from footer: '+str(mv.exptime)+' s')
    #now display the first frame
    processframe()
    displayFrame(autoscale=True,markstars=False)
    #Load calibration frames and set up
    processLog.log("Please load dark, flat, and dark for flat files",1)
    #Select stars:
    selectstars()
    #mv.spe.close() #In real version we'll close spe
    win.setupApsizeMenu()
        

#Determine the exposuretime of a SPE file without a footer
def getexptime(thisspe):
    #Input open SPE file
    #don't read lots of frames in large files
    numtoread = min([thisspe.get_num_frames(),11])
    tstamps = np.zeros(numtoread)
    for f in range(numtoread): 
        tstamps[f] = mv.spe.get_frame(f)[1]['time_stamp_exposure_started']
    timediff = tstamps[1:numtoread]-tstamps[:numtoread-1]
    return np.round(np.median(timediff/1e6))


#Define all the stuff that needs to be done to each incoming frame
def processframe(i=0):
    (thisframe,thistime) = mv.spe.get_frame(i)
    #calibrate (doesn't do anything if calibration frames are not available):
    if mv.dark != []: thisframe=(thisframe-mv.dark)
    if mv.flat != []: thisframe=thisframe/mv.flat
    #read in frame
    mv.img=np.transpose(thisframe)
    backgroundmed,backgroundvar=charbackground()

    #Replace if this frame already exists, otherwise append
    if i <= mv.framenum: #replace
        #log('Re-processing frame '+str(i)+' of '+str(mv.framenum))
        mv.rawtimes[i]=thistime['time_stamp_exposure_started']
        mv.backmed[i]=backgroundmed
        mv.backvar[i]=backgroundvar
    else: #append
        #log('Processing frame '+str(i)+' of '+str(mv.framenum))
        mv.rawtimes.append(thistime['time_stamp_exposure_started'])
        mv.backmed.append(backgroundmed)
        mv.backvar.append(backgroundvar)
    
    mv.framenum=i
    
#Function to characterize the background to find stellar centroids accurately
#This should be done for each frame as it's read in
def charbackground():
    """Characterize the image background, median and variance

    for frame currenly held in img  
    """
    backgroundmed = biweight_location(mv.img)
    backgroundvar = biweight_midvariance(mv.img)
    return backgroundmed, backgroundvar        
 
#show the image to the widget
def displayFrame(autoscale=False,markstars=True):
    """Display an RBG image
    
    i is index to display
    Autoscale optional.
    Return nothing.
    """
    #Make sure i is in range
    '''
    if autoscale:
        #lowlevel=np.min(thisimg[thisimg > 0])
        lowlevel=np.min(mv.displayimg)
        if np.sum(mv.displayimg > 0) > 100:
            lowlevel=np.percentile(mv.displayimg[mv.displayimg > 0],3)
        highlevel=np.max(mv.displayimg)-1
        w5.setImage(np.array(mv.displayimg),autoRange=True,levels=[lowlevel,highlevel],)
    else:
        w5.setImage(np.array(mv.displayimg),autoRange=False,autoLevels=False)
    #Draw position circles:
    if markstars and len(mv.stars) > 0:
        targs.setData([p[0] for p in mv.stars[mv.framenum]],[p[1] for p in mv.stars[mv.framenum]])
        targs.setSize(2.*mv.apsizes[mv.apsizeindex])
        targs.setPen(pencolors[0:mv.numstars])
    '''
    imageDisplay.displayImage(mv.img)
    if markstars and len(mv.stars) > 0:
        imageDisplay.displayApertures(mv.stars[mv.framenum],mv.apsizes[mv.apsizeindex])
        
        
def selectstars():
    '''Select stars in the current frame.
    
    Click to select any number in the first image.
    Click to select numstars in later images to get following back on track.
    '''
    mv.selectingstars = True
    
def gaussian(x, A, sigma):
    #Define a gaussian for finding FWHM
    return A*np.exp(-(x)**2/(2.*sigma**2))

def improvecoords(x,y,i=mv.framenum,pixdist=10,fwhm=4.0,sigma=5.):
    """Improve stellar centroid position from guess value. (one at a time)

    #return the adjustment than needs to be made in x and y directions
    #also calculate the FWHM seeing
    """
    #x=(1024/binning)-x
    #y=(1024/binning)-y
    #Keep track of motion
    delta = np.zeros(2)
    #Get image subregion around guess position
    #Need to be careful not to ask for out-of-range indexes near a border
    x0=x-pixdist
    y0=y-pixdist
    xdist=2*pixdist
    ydist=2*pixdist
    if x0 < 0: #if near the left edge
        x0 = 0 #subregion from near given position
        delta[0] += pixdist-x #adjust delta accordingly
    if y0 < 0: #same in the y direction
        y0 = 0
        delta[1] += pixdist-y
    if x+pixdist > mv.img.shape[0]:
        xdist = mv.img.shape[0]-x+pixdist
    if y+pixdist > mv.img.shape[1]:
        ydist = mv.img.shape[1]-y+pixdist
    subdata=mv.img[x0:x0+xdist,y0:y0+ydist]
    #print subdata.shape
    sources = daofind(subdata - mv.backmed[i], sigma*mv.backvar[i], fwhm,
                      sharplo=0.1, sharphi=1.5, roundlo=-2.0, roundhi=2.0)
    #From what I can tell, daofind returns x and y swapped, so fix it
    returnedx = sources['ycentroid']
    returnedy = sources['xcentroid']
    
    thisseeing = np.nan 
    
    if len(sources) != 0:
        strongsignal= np.argmax(sources['peak'])
        delta[0]+=returnedx[strongsignal]-pixdist
        delta[1]+=returnedy[strongsignal]-pixdist
        #Fit with a gaussian
        seeingdata = subdata.flatten() - mv.backmed[i]
        dist = []
        for j in np.arange(subdata.shape[1])+0.5:
            for k in np.arange(subdata.shape[0])+0.5:
                dist.append(np.sqrt((returnedy[strongsignal]-k)**2.
                            +(returnedx[strongsignal]-j)**2.))
        dist=np.array(dist).flatten()#distance between new coord and pixel centers
        #plt.scatter(dist,seeingdata)
        try: #ignores error if max iterations is hit        
            p0=[1000.,4.]#initial guesses
            popt,_  = curve_fit(gaussian,np.append(dist,dist*-1.),np.append(seeingdata,seeingdata),p0=p0)
            thisseeing = np.abs(popt[-1])*2.3548
            #plt.plot(np.arange(0,10,.1),gaussian(np.arange(0,10,.1),popt[0],popt[1]))
        except RuntimeError:
            print "ERROR: gaussian fit did not converge for a star in frame "+str(i)
        #plt.show()
    else:
        delta=np.zeros(2)

    #also measure the seeing in this step:


    #check that unique source found
    '''
    if len(sources) == 0:
        processLog.log("Frame #"+str(i),1)
        processLog.log("WARNING: no sources found in searched region near ({:.2f}, {:.2f}).".format(x,y))
        #delta = [0,0] in this case
    else:
        if len(sources) > 1:
            processLog.log("Frame #"+str(i),1)
            processLog.log("WARNING: non-unique solution found for target near ({:.2f}, {:.2f}).".format(x,y))
            processLog.log(str(len(sources))+" signals in window.  Using brightest.")
        #Take brightest star found
    '''

    #handle stars that were not found #Move this outside this function
    """
    if [0,0] in delta and follow:
        meandeltax=np.mean(delta[np.where(delta[:,0] != 0),0])
        meandeltay=np.mean(delta[np.where(delta[:,1] != 0),1])
        delta[np.where(delta[:,0] == 0)] += [meandeltax,meandeltay]
    """

    return delta[0],delta[1],thisseeing
















#### STAGE 2 ####

#Aperture details (provide a way to change these!)
r_in = 16.  #inner sky annulus radius #change in terms of binning eventually
r_out = 24. #outer sky annulus radius #change in terms of binning eventually
def setApSize(size):
    processLog.log("Aperture size set to "+str(size)+" pixels.",1)
    #processLog.log("(Updates on next frame.)")
    if size in mv.apsizes:
        mv.apsizeindex=np.where(mv.apsizes == size)[0][0]
        if markstars and len(mv.stars) > 0:
            imageDisplay.displayApertures(mv.stars[mv.framenum],size)
        if mv.stage > 1: 
            updatelcs(i=mv.framenum)


def setCompStar(s):
    mv.compstar = s
    processLog.log("Now dividing by comparsion star #"+str(s),1)
    updatelcs(mv.framenum)





#Run the stage 2 loop
def stage2():
    stagechange(2)
    #Add plot items for raw counts panel to plot
    for splot in rawcounts: w7.addItem(splot)
    for splot in seeingplots: w9.addItem(splot)
    #Make stars array an array of arrays of star coord arrays (yikes)
    # i.e, it needs to get pushed a level deeper
    mv.stars=[mv.stars]
    #same with seeing
    mv.seeing=np.array([mv.seeing])
    #Run photometry on the first frame
    dophot(0)
    updatelcs(i=0)
    updatehack()
    #Start timer that looks for new data
    timer2.start(min(mv.exptime*1000.,5000.))# shorter of exptime and 5 sec
    timer3.start(1.*60*1000)#update every 1 minutes
    #This currently freezes up the UI.  Need to thread, but not enough time
    #to implement this currently.  Use a hack for now
    '''
    #Run the loop:
    mv.fsize_spe_old = 0
    while not mv.hasFooter:
        #Update only if there's new data
        fsize_spe_new = os.path.getsize(mv.spefile)
        if fsize_spe_new > mv.fsize_spe_old:
            mv.spe = read_spe.File(mv.spefile)
            mv.numframes = spe.get_num_frames()
            processLog.log('Processing frames '+str(mv.framenum)+'-'+str(mv.numframes),1)
            while mv.framenum < mv.numframes:
                nextframe()
            if hasattr(mv.spe, 'footer_metadata'): 
                mv.hasFooter = True
                processLog.log('SPE footer detected. Data acquisition complete.',2)
                stagechange(3)
            mv.spe.close()
        mv.fsize_spe_old = fsize_spe_new
    '''


def updatehack():
    #Only look for new data if not currently processing new data
    if not timer.isActive():
        #Update only if there's new data
        fsize_spe_new = os.path.getsize(mv.spefile)
        
        if fsize_spe_new > mv.fsize_spe_old and mv.stage ==2:
            mv.spe = read_spe.File(mv.spefile)
            mv.numframes = mv.spe.get_num_frames()
            if mv.framenum+1==mv.numframes-1:processLog.log('Processing frame '+str(mv.framenum+1))
            else: processLog.log('Processing frames '+str(mv.framenum+1)+'-'+str(mv.numframes-1),1)        
            timer.start(100)
            #Update plots
            updatelcs(i=mv.framenum)
            if hasattr(mv.spe, 'footer_metadata'): 
                mv.hasFooter = True
                timer3.stop()
        mv.fsize_spe_old = fsize_spe_new

def nextframehack():
    #call nextframe until you're caught up
    nextframe()
    updatelcs(i=mv.framenum)
    if mv.framenum >= mv.numframes-1:
        timer.stop()
        updateft(i=mv.framenum)
        if mv.hasFooter:
            processLog.log('SPE footer detected. Data acquisition complete.',2)
            stagechange(3)
            processLog.log("Image processing complete",2)
            writetimestamps()
            displayFrame(autoscale=True)
            
        mv.spe.close()

#This timer catches up on photometry
timer = pg.QtCore.QTimer()#set up timer to avoid while loop
timer.timeout.connect(nextframehack)

#This timer checks for new data
timer2 = pg.QtCore.QTimer()
timer2.timeout.connect(updatehack)



#For demo purposes, read in the next frame of the spe file each time this is called
def nextframe():
    #if mv.stage == 2:
    oldcoords = mv.stars[mv.framenum]
    processframe(i=mv.framenum+1) #Frame num increases here.
    newcoords=[]
    newseeing=[]
    for coord in oldcoords:
        dx,dy,thisseeing = improvecoords(coord[0],coord[1],i=mv.framenum)
        newcoords.append([np.floor(coord[0])+.5+dx,np.floor(coord[1])+.5+dy]) 
        newseeing.append(thisseeing)
    mv.stars.append(newcoords)
    mv.seeing = np.append(mv.seeing,[newseeing],axis=0)
    #Show the frame
    displayFrame(autoscale=True,markstars=True)
    #Perform photometry
    dophot(i=mv.framenum)
    #Update light curves
    #updatelcs(i=mv.framenum) #only after all the new photometry is done.
    


    

def dophot(i):
    '''Do photometric measurements.
    
    Stars have been selected.  Do aperture photometry on given frame
    '''
    #print "dophot(i) called with i="+str(i)
    #Do the aperture photometry

    #The aperture_photometry() function can do many stars at once
    #But you must first do a background subtraction


    #We're going to save a lot of information in this step:
    #Total counts and uncertainty for every aperture size for every star
    #And eventually for every frame...

    #Note that the photometry package seems to reference x and y coords
    #as the tranpose of what we've been using.  Switch the order here:
    coords = [star[::-1] for star in mv.stars[i]]
    thisphotometry = np.zeros((len(coords),len(mv.apsizes)))
    for n in range(mv.numstars):
        #Loop through the stars in the image
        #annulus_aperture = CircularAnnulus(coords[n], r_in=r_in, r_out=r_out)
        #print aperture_photometry(mv.img[i],annulus_aperture).keys()
        #background_mean = aperture_photometry(mv.img[i],annulus_aperture)['aperture_sum'][0]/annulus_aperture.area() 
        #NOTE on the above line: This should really be a median!
        #Issue 161 on photutils https://github.com/astropy/photutils/issues/161 is open as of 09/28/15
        gain = 12.63 #? From PI Certificate of Performance for "traditional 5MHz gain."  Confirm this value!
        #loop through aperture sizes
        for j,size in enumerate(mv.apsizes):
            aperture = CircularAperture(np.array(coords[n]), r=size) 
            #phot = aperture_photometry(x-background_mean,aperture,error=backgroundvar,gain=gain)
            #Why am I getting negative numbers?
            #phot = aperture_photometry(mv.img[i]-np.median(mv.img),aperture)
            phot = aperture_photometry(mv.img-mv.backmed[i],aperture)
            thisphotometry[n,j]=phot['aperture_sum'][0]
    #print "photometry ",thisphotometry
    if i == 0:
        mv.photresults = np.array([thisphotometry])
    else:
        mv.photresults = np.append(mv.photresults,[thisphotometry],axis=0)
    #print "photresults dimensions are "+str(mv.photresults.shape)
    #yay.  This deserves to all be checked very carefully, especially since the gain only affects uncertainty and not overall counts.


#Allow different kernel types:
kerneltypes = ['Uniform','Epanechnikov']

#set up a dialog to change the kernel details:
class KernelDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(KernelDialog, self).__init__(parent)
        typeLabel = QtGui.QLabel("Kernel &type")
        self.typeEdit = QtGui.QComboBox()
        self.typeEdit.addItems(kerneltypes)
        #self.typeEdit.setCurrentIndex(currentind)
        typeLabel.setBuddy(self.typeEdit)
        widthLabel = QtGui.QLabel("Kernel &width")
        self.widthEdit = QtGui.QSpinBox()
        self.widthEdit.setMinimum(2)
        self.widthEdit.setMaximum(200)
        #self.widthEdit.setValue(currentwidth)
        widthLabel.setBuddy(self.widthEdit)
        self.buttons = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
                                        QtCore.Qt.Horizontal, self)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        
        grid = QtGui.QGridLayout(self)
        grid.addWidget(typeLabel,0,0)
        grid.addWidget(self.typeEdit,0,1)
        grid.addWidget(widthLabel,1,0)
        grid.addWidget(self.widthEdit,1,1)
        grid.addWidget(self.buttons, 3, 0)
        self.setLayout(grid)

        self.setWindowTitle("Define Smoothing Kernel")
    def kernelFormat(self):
        kerneltype=int(self.typeEdit.currentIndex())
        width=int(self.widthEdit.value())
        return (kerneltype,width)
        
    @staticmethod
    def getKernelFormat(parent = None):
        dialog = KernelDialog(parent)
        result = dialog.exec_()
        kerneltype,width = dialog.kernelFormat()
        return (kerneltype,width, result == QtGui.QDialog.Accepted)


#set up a class that holds all the smoothing kernel information
class smoothingkernel:
    """Holds all smoothing kernel info"""
    kerneltype = 0
    width = 10 #points
    kernel=[]
    types = kerneltypes
    def setkernel(self,kerneltype,width):
        if kerneltype == 1: #Epanechnikov
            u=(2.*np.arange(width)/(float(width)-1.))-0.5
            self.kernel = 0.75*(1.-u**2.)
            self.kernel /= np.sum(self.kernel)
            processLog.log("Using "+self.types[kerneltype]+" smoothing kernel of width "+str(width))
        elif kerneltype == 0: #Uniform
            self.kernel = np.ones(width)/float(width)
            processLog.log("Using "+self.types[kerneltype]+" smoothing kernel of width "+str(width))
    def openKernelDialog(self):
            dkerneltype,dwidth,daccepted = KernelDialog.getKernelFormat()
            if daccepted and (dkerneltype in range(len(kerneltypes))) and (dwidth > 1): 
                self.setkernel(dkerneltype,dwidth)
    def __init__(self):
        self.setkernel(0,10)

#set up the kernel object
kernel=smoothingkernel()


#Update display.
def updatelcs(i):
    #Identify which points to include/exclude, up to frame i
    goodmask=np.ones(i+1, np.bool)
    goodmask[mv.bad] = False
    badmask = np.zeros(i+1, np.bool)
    badmask[mv.bad] = True
    targdivided = mv.photresults[:i+1,0,mv.apsizeindex]/mv.photresults[:i+1,mv.compstar,mv.apsizeindex]
    times = np.arange(i+1)#Multiply by exptime for timestamps
    goodfluxnorm=targdivided[goodmask[:i+1]]/np.abs(np.mean(targdivided[goodmask[:i+1]]))
    s1.setData(mv.exptime*times[goodmask[:i+1]],goodfluxnorm)
    #s2.setData(times[badmask[:i]],targdivided[badmask[:i]])
    l1.setData(mv.exptime*times[goodmask[:i+1]],goodfluxnorm)
    #sl1.setData(times[goodmask[:i]],fluxsmoothed[goodmask[:i]])
    #Raw Counts:
    for j,splot in enumerate(rawcounts): splot.setData(mv.exptime*times,mv.photresults[:,j,mv.apsizeindex])
    #Seeing:
    for j,splot in enumerate(seeingplots[::-1]): splot.setData(mv.exptime*times,mv.seeing[:,j])
    #Sky brightness
    sky.setData(mv.exptime*times,mv.backmed)

def updateftfromtimer():
    updateft(i=mv.framenum)

def updateft(i=mv.framenum):
        goodmask=np.ones(i+1, np.bool)
        goodmask[mv.bad] = False
        targdivided = mv.photresults[:i+1,0,mv.apsizeindex]/mv.photresults[:i+1,mv.compstar,mv.apsizeindex]
        goodfluxnorm=targdivided[goodmask[:i+1]]/np.abs(np.mean(targdivided[goodmask[:i+1]]))
        times = np.arange(i+1)#Multiply by exptime for timestamps
        #Fourier Transform   and smoothed lc  
        ftplot.calcft(times*mv.exptime,goodfluxnorm-1.,mv.exptime)
        
        #Smoothed LC
        if goodmask.sum() > 2:
            #This all requires at least two points
            #Only update once per file read-in
            interped = interp1d(mv.exptime*times[goodmask[:i+1]],goodfluxnorm-1.)
            xnew = np.arange(mv.exptime*min(times[goodmask[:i]]),mv.exptime*max(times[goodmask[:i+1]]),mv.exptime)
            ynew = interped(xnew)

            #Smoothed LC
            #Update if there are enough points:
            if len(ynew) > kernel.width:
                fluxsmoothed=np.convolve(ynew,kernel.kernel,mode='same')
                ss1.setData(xnew,fluxsmoothed)


#This timer recomputes the FT and smoothed lc infrequently
timer3 = pg.QtCore.QTimer()
timer3.timeout.connect(updateftfromtimer)

''' Not implemented yet!
#To keep the GUI from locking up, computationally intensive processes must
#be done in a thread.  Set up that thread here:
class Stage2Thread(QtCore.QThread):

    setTime = QtCore.pyqtSignal(int,int)
    iteration = QtCore.pyqtSignal(threading.Event, int)

    def run(self):

        self.setTime.emit(0,300)
        for i in range(300):
            time.sleep(0.05)
            event = threading.Event()
            self.iteration.emit(event, i)
            event.wait()

'''


#Write timestamps
def writetimestamps():
    fpath_csv = os.path.splitext(mv.spefile)[0]+'_timestamps.csv'
    processLog.log("Writing absolute timestamps to file "+fpath_csv,2)
    if hasattr(mv.spe, 'footer_metadata'):
        footer_metadata = BeautifulSoup(mv.spe.footer_metadata, "xml")
        trigger_response = footer_metadata.find(name='TriggerResponse').text
        ts_begin = footer_metadata.find(name='TimeStamp', event='ExposureStarted').attrs['absoluteTime']
        dt_begin = dateutil.parser.parse(ts_begin)
        ticks_per_second = int(footer_metadata.find(name='TimeStamp', event='ExposureStarted').attrs['resolution'])
    else:
        processLog.log(("No XML footer metadata.\n" +
               "Unknown trigger response.\n" +
               "Using file creation time as absolute timestamp.\n" +
               "Assuming 1E6 ticks per seconds."),3)
        trigger_response = ""
        dt_begin = dt.datetime.utcfromtimestamp(os.path.getctime(fpath_spe))
        ticks_per_second = 1E6
    idx_metadata_map = {}
    for idx in xrange(mv.spe.get_num_frames()):
        (frame, metadata) = mv.spe.get_frame(idx)
        idx_metadata_map[idx] = metadata
    df_metadata = pd.DataFrame.from_dict(idx_metadata_map, orient='index')
    df_metadata = df_metadata.set_index(keys='frame_tracking_number')
    df_metadata = df_metadata[['time_stamp_exposure_started', 'time_stamp_exposure_ended']].applymap(lambda x: x / ticks_per_second)
    df_metadata = df_metadata[['time_stamp_exposure_started', 'time_stamp_exposure_ended']].applymap(lambda x : dt_begin + dt.timedelta(seconds=x))
    df_metadata[['diff_time_stamp_exposure_started', 'diff_time_stamp_exposure_ended']] = df_metadata - df_metadata.shift()
    processLog.log("Trigger response = {tr}".format(tr=trigger_response))
    processLog.log("Absolute timestamp = {dt_begin}".format(dt_begin=dt_begin))
    processLog.log("Ticks per second = {tps}".format(tps=ticks_per_second))
    df_metadata.head()
    
    # Write out as CSV to source directory of SPE file.
    df_metadata.to_csv(fpath_csv, quoting=csv.QUOTE_NONNUMERIC)
    saveScreenshot()






## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        if len(sys.argv) > 1:
            mv.defaultdir = sys.argv[1]
        QtGui.QApplication.instance().exec_()
