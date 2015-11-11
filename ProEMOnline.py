# -*- coding: utf-8 -*-
"""
This scripts sets an initial layout for the ProEMOnline software.  It uses the
PyQtGraph dockarea system and was designed from the dockarea.py example.


Keaton wrote this.
"""


#Import everything you'll need
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pickle #for saving layouts
from glob import glob
from scipy import stats
from scipy.fftpack import fft,fftfreq
import os
import csv
import sys
import time
import datetime as dt
import dateutil.parser
from astropy.io import fits
from astroML.fourier import FT_continuous
from scipy.interpolate import interp1d
import scipy.ndimage.filters as filters
from astropy.stats import biweight_location, biweight_midvariance
from photutils import daofind
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from pyqtgraph.dockarea import *
from bs4 import BeautifulSoup
# Local modules.
import read_spe



#### BEGIN PROGRAM ####


#The organization and behavoir of the program are as follows:
#This program operates in four stages.
#Stage 0 - Program Initialized, waiting to open SPE file.
#Stage 1 - SPE file open, stars are being selected
#Stage 2 - Online data reduction and aperture photometry/plotting is being done.
#Stage 3 - End of data acquisition detected. Final data written to file.  Timestamps verified.  Log saved.  Weather/time log data saved.
# -> revert back to Stage 0.
stage=0 #start at 0
def stagechange(num):
    global stage
    if num in range(4):
        log("Program stage = "+str(num),1)
        stage=num
    else: log("Attempt to change stage to invalid value ("+str(num)+")",3)



#### STAGE 0 ####
#Set up the general GUI aspects
defaultdir = 'D:/sync_to_White_Dwarf_Archive/'#where to search for SPE files


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
        
        #Run Photometry
        updateFT = QtGui.QAction('&Update FT', self)
        updateFT.setShortcut('Ctrl+U')
        updateFT.setStatusTip('Update Fourier Transform with Current Light Curve')
        updateFT.triggered.connect(self.updateFTfunct)
        
        #Load dark for science frames
        loadDark = QtGui.QAction('Load Dark SPE', self)
        loadDark.setStatusTip('Open SPE Calibrations for Dark Subtracting Science Images')
        loadDark.triggered.connect(self.openDark)
        
        #Load dark for flat frames
        loadDarkForFlats = QtGui.QAction('Load Dark for Flats', self)
        loadDarkForFlats.setStatusTip('Open SPE Calibrations for Dark Subtracting Flat Images')     
        loadDarkForFlats.triggered.connect(self.openDarkForFlats)
        
        #Load flat
        loadFlat = QtGui.QAction('Load Flat SPE', self)
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
        
        #Menubar
        menubar = self.menuBar()
        #File Menu
        fileMenu = menubar.addMenu('File')
        fileMenu.addAction(openFile)
        fileMenu.addAction(runPhot)
        fileMenu.addAction(updateFT)
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
        #Layout Menu
        layoutMenu = menubar.addMenu('Layout')
        layoutMenu.addAction(saveLayout)
        layoutMenu.addAction(loadLayout)
        

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
        global defaultdir,rundir
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open SPE file', 
                defaultdir,filter='Data (*.spe)'))
        if fname[-4:]=='.spe':
            log("Opening file "+fname,1)
            #Set the default directory to a couple levels up from this file
            rundir = os.path.dirname(fname)
            defaultdir = os.path.dirname(rundir)
            #set target log text as filename to start
            targetEdit.setText(os.path.basename(fname)[:-4])
            
            #This needs to trigger a major chain of events
            stage1(fname)
        else: log("Invalid file type (must be SPE).",3)
        
    #Update the FT at user's command
    def updateFTfunct(self):
        global framenum
        updateft(i=framenum)
        
    #Load Dark frames
    def openDark(self):
        global dark, darkExists
        fname = str(QtGui.QFileDialog.getOpenFileName(self, 'Open dark file', 
                defaultdir,filter='Data (*.spe *.fits)'))
        if fname[-4:]=='.spe':
            log("Opening dark file "+fname,1)
            dspe = read_spe.File(fname)
            num_darks=dspe.get_num_frames()
            #get all frames in SPE file
            #stack as 3D numpy array
            (frames,_)=dspe.get_frame(0)
            frames=np.array([frames])
            for i in range(1,num_darks):
                (thisframe,_)=dspe.get_frame(i)
                frames=np.concatenate((frames,[thisframe]),0)
            dark=np.median(frames,axis=0)
            darkExists = True
            log("Mean dark counts: "+str(np.mean(dark)))
            processframe()
            displayFrame(autoscale=True,markstars=False)
            
            #Write out master dark file as fits     
            #Set up header
            prihdr = fits.Header()
            prihdr['OBJECT'] = 'dark'
            prihdr['IMAGETYP'] = 'dark'
            prihdr['REDUCED'] = dt.datetime.now().isoformat()
            prihdr['COMMENT'] = "Reduced by Keaton Bell's OLD MAID Software"
            if hasattr(dspe, 'footer_metadata'):
                footer_metadata = BeautifulSoup(dspe.footer_metadata, "xml")
                ts_begin = footer_metadata.find(name='TimeStamp', event='ExposureStarted').attrs['absoluteTime']
                dt_begin = dateutil.parser.parse(ts_begin)
                prihdr['TICKRATE'] = int(footer_metadata.find(name='TimeStamp', event='ExposureStarted').attrs['resolution'])
                prihdr['DATE-OBS'] = str(dt_begin.isoformat())
                prihdr['XBINNING'] = footer_metadata.find(name="SensorMapping").attrs['xBinning']
                prihdr['YBINNING'] = footer_metadata.find(name="SensorMapping").attrs['yBinning']
                prihdr['INSTRUME'] = footer_metadata.find(name="Camera").attrs['model']
                prihdr['TRIGGER'] = footer_metadata.find(name='TriggerResponse').text
                prihdr['COMMENT'] = "SPE file has footer metadata"
                prihdr['EXPTIME'] = str(float(footer_metadata.find(name='ExposureTime').text)/1000.)
                #prihdr['SOFTWARE'] = footer_metadata.find(name='Origin')
                prihdr['SHUTTER'] = footer_metadata.find(name='Mode').text
                if footer_metadata.find(name='Mode').text != 'AlwaysClosed':
                    prihdr['WARNING'] = 'Shutter not closed for dark frame.'
                    log("Shutter not closed for dark frame.",3)
            else:
                prihdr['WARNING'] = "No XML footer metadata."
                log("No XML footer metadata.",3)
            #Set up fits object
            hdu = fits.PrimaryHDU(dark,header=prihdr)
            darkpath = os.path.dirname(fname)
            fitsfilename = 'master_'+os.path.basename(fname).split('.spe')[0]+'.fits'
            log("Writing master dark as "+fitsfilename)
            hdu.writeto(os.path.join(darkpath, fitsfilename),clobber=True)
            #Close SPE
            dspe.close()
            
        elif fname[-5:]=='.fits':
            log("Opening dark file "+fname,1)
            hdulist = fits.open(fname)
            prihdr = hdulist[0].header
            dark=hdulist[0].data
            log("Mean dark counts: "+str(np.mean(dark)))
            darkExists = True
            processframe()
            displayFrame(autoscale=True,markstars=False)
            hdulist.close()
        else: log("Invalid file type (must be SPE).",3)
        
        
    #Load Dark frames for flat calibration
    def openDarkForFlats(self):
        global darkForFlat, darkExistsForFlat
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open SPE dark for flat calibration', 
                defaultdir,filter='Data (*.spe)')
        if fname[-4:]=='.spe':
            log("Opening dark file "+fname+" for flat calibration.",1)
            dspe = read_spe.File(str(fname))
            num_darks=dspe.get_num_frames()
            #get all frames in SPE file
            #stack as 3D numpy array
            (frames,_)=dspe.get_frame(0)
            frames=np.array([frames])
            for i in range(1,num_darks):
                (thisframe,_)=dspe.get_frame(i)
                frames=np.concatenate((frames,[thisframe]),0)
            darkForFlat=np.median(frames,axis=0)
            darkExistsForFlat = True
            log("Mean dark counts for flat: "+str(np.mean(darkForFlat)))
            dspe.close()
        else: log("Invalid file type (must be SPE).",3)        
        
        
    #Load Flat frames
    def openFlat(self):
        global flat, flatExists
        if darkExistsForFlat == False:
            log("Import dark for reducting flats before importing flat file.",3)
        else:
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Open SPE flat file', 
                    defaultdir,filter='Data (*.spe)')
            if fname[-4:]=='.spe':
                log("Opening flat file "+fname,1)
                fspe = read_spe.File(str(fname))
                num_flats=fspe.get_num_frames()
                #get all frames in SPE file
                #stack as 3D numpy array
                (frames,_)=fspe.get_frame(0)
                modes=[]
                frames = frames - darkForFlat
                modes.append(stats.mode(frames.flatten())[0][0])
                frames=np.array([frames/modes[0]])
                for i in range(1,num_flats):
                    (thisframe,_)=fspe.get_frame(i)
                    thisframe = thisframe-darkForFlat
                    modes.append(stats.mode(thisframe.flatten())[0][0])
                    frames=np.concatenate((frames,[thisframe/modes[i]]),0)
                flat=np.median(frames,axis=0)
                flatExists=True
                log("Median flat counts: "+str(np.median(modes)))
                fspe.close()
                processframe()
                displayFrame(autoscale=True,markstars=False)
            else: log("Invalid file type (must be SPE).",3)
    
    
    #Restore previously "bad" points
    def restorePts(self):
        global bad
        log("Deselecting "+str(len(bad))+" points.")
        bad=[]
        updatelcs(i=framenum)
    
    #Set up aperture size menu options
    def setupApsizeMenu(self): 
        for size in apsizes:
            self.changeApertureMenu.addAction(str(size)+' pixels',lambda s=size: setApSize(s))
    
    #Set up comp star selection menu options
    def addCompStarOption(self,i): 
        self.changeCompStarMenu.addAction('Comp Star #'+str(i),lambda s=i: setCompStar(s))
            
    
    #Run Photometry
    def run(self):
        #Do aperture photometry on selected stars
        global numstars, selectingstars
        if stage == 1:
            if len(stars) == 0:
                log("No stars selected.  Select stars before running.",3)
            else:
                numstars = len(stars)
                #Write original coordinates to phot_coords.orig
                f = open(rundir+'/phot_coords.orig', 'w')
                for star in stars:
                    f.write('{:.2f} {:.2f}\n'.format(star[0],star[1]))
                f.close()
                selectingstars=False
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
processLog = QtGui.QTextEdit()
processLog.setReadOnly(True)
w2.addWidget(processLog, 0, 0, 6, 1)
d2.addWidget(w2)
# This widget need special functions to get messages:
def log(text,level=0):
    '''log messages to the process log and log file
    
    text is the message for the log
    level indicated how important it is:
    level=0: Routine background process: gray text;
    level=1: User influenced action: black text;
    level=2: Major change: bold black;
    level=3: Warning message: bold red; 
    '''
    text=str(text)
    colors = ['darkgray','black','black','red']
    prefix = ['','','','WARNING: ']
    fontweight = [50,50,75,75]
    if level in range(4):
        processLog.setTextColor(QtGui.QColor(colors[level]))
        processLog.setFontWeight(fontweight[level])
        processLog.append(prefix[level]+text)
    else: log('Level assigned to message "'+text+'" out of range.',level=3)


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
    global bad
    for p in points:
        if p.pos()[0]/exptime in bad:
            bad.remove(p.pos()[0]/exptime)
        else:
            bad.append(p.pos()[0]/exptime)
    updatelcs(i=framenum)
    
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
w9 = pg.PlotWidget(title="Seeing (not yet implemented)",labels={'left': 'seeing (")', 'bottom': 'time (s)'})
seeing = pg.PlotCurveItem()
w9.addItem(seeing) #Placeholder for now
d9.addWidget(w9)


## Fourier Transform
w3 = pg.PlotWidget(title="Fourier Transform",labels={'left': 'amplitude (mma)', 'bottom': 'freq (muHz)'})
ft = w3.plot(pen='y')
d3.addWidget(w3)


## Image
w5 = pg.ImageView()
w5.ui.roiBtn.hide()
#w5.ui.normBtn.hide() #Causing trouble on windows
#Define function for selecting stars. (must be defined before linking the click action)
def click(event):#Linked to image click event
    global stars
    if event.button() == 1 and selectingstars:
        event.accept()
        pos = event.pos()
        #x and y are swapped in the GUI!
        x=pos.x()
        y=pos.y()
        #log('Clicked at ({:.2f}, {:.2f})'.format(x,y),level=0)
        #improve coordinates
        dx,dy = improvecoords(x,y)
        #round originals so original position *within* pixel doesn't affect answer
        newcoords=[np.floor(x)+dx,np.floor(y)+dy]
        stars.append(newcoords)
        #make menuoption for comp star selection
        if len(stars) > 1: win.addCompStarOption(len(stars)-1)
        #Mark stars in image display
        targs.setData([p[0] for p in stars],[p[1] for p in stars])
        targs.setPen(pencolors[0:len(stars)])
        #Set up plot for raw counts in second panel:
        rawcounts.append(pg.ScatterPlotItem(pen=pencolors[len(stars)-1],symbol='o',size=1))
        log('Star selected at ({:.2f}, {:.2f})'.format(newcoords[0],newcoords[1]),level=2)
        
    elif event.button() == 2: 
        event.accept()#Passed on to other functionality if not accepted.
        print "RIGHT!"


w5.getImageItem().mouseClickEvent = click #Function defined below
#w5.keyPressEvent = moveCircles # Seems to be the right thing for detecting frame changes,
#But I can't connect to it without overriding other behavior.  May need to subclass this.


#Set up plot for apertures around stars
#print QtGui.QColor.colorNames() for available names.
stringcolors=['red','green','blue','magenta','orange','yellow',
              'darkred','darkgreen','darkblue','darkmagenta','darkorange','darkgoldenrod',
              'hotpink','seagreen','skyblue','salmon','brown','lightyellow']
pencolors = [pg.mkPen(QtGui.QColor(c), width=3) for c in stringcolors]
targs = pg.ScatterPlotItem(brush=None, pen=pencolors[0],symbol='o',pxMode=False,size=8)
w5.addItem(targs)
#Add widget to dock

d5.addWidget(w5)


## Show the program!
win.show()
win.raise_()

#win.activateWindow()













# I think everything is set up enough to start doing stuff
# Send initial message to process log.
log("ProEMOnline initialized",2)
#log("Development version.  Do not trust.",3)
stagechange(0)
log("Open SPE file to begin analysis.",1)



#### STAGE 1 ####

# Stage 1 starts when a SPE file is loaded.
# It's the "getting everything set up" stage
# Since the SPE file is loaded by the menu action, this will be one big
# function that is called on the new image.

#First define all the variables everything will need access to:
#These will be called into action as global variables.

#SPE Filename
spefile = ''
#SPE Data
spe=[]
#SPE file directory
rundir=''
#Does SPE have a footer?
hasFooter=False
#Number of frames in currently read spe file
numframes=0
#Exposure time for science frames
exptime=1. #If it can't be figured out, plots are in terms of frame #
#Dark data
dark = []
darkExists=False
darkForFlat = []
darkExistsForFlat=False
#Flat data
flat = []
flatExists=False
#Number of last *reduced* (photometry measures) frame
framenum=-1 #none yet
#Flag to indicate whether we are currently selecting stars in the frame:
selectingstars = False
#Number of stars to do photometry on (target first)
numstars = 0 #0 means we haven't selected stars yet.
#Star coords
stars = [] #list of list of list of coords
#Image data:
img=[] #only hold current image to save tiem
#And another version to look nice
displayimg=[] #only hold current image to save tiem
#Keep track of "Bad" points
bad=[]
#Elapsed timestamps
rawtimes=[] #start of timestamp
#Search radius (box for now), improve later
pixdist=10 #(super)pixels
#List of median background counts:
backmed=[]
#List of background variances
backvar=[]
#Binning
binning=4

def stage1(fname):
    #Load SPE File

    #Access needed global vars
    global spefile,spe,binning,exptime,dark,flat
    #Announce Stage 1    
    stagechange(1)
    #Record SPE filename this once
    spefile = fname
    #Read in SPE data
    spe = read_spe.File(spefile)
    binning = 1024/spe.get_frame(0)[0].shape[0]
    log(str(spe.get_num_frames()) + ' frames read in.')
    exptime=getexptime(spe)
    log('Inferred exposure time: '+str(exptime)+' s')
    if hasattr(spe, 'footer_metadata'): 
        #log('SPE file has footer.')
        exptime=np.round(float(BeautifulSoup(spe.footer_metadata, "xml").find(name='ExposureTime').text)/1000.)
        #log('Exposute time from footer: '+str(exptime)+' s')
    #now display the first frame
    processframe()
    displayFrame(autoscale=True,markstars=False)
    #Load calibration frames and set up
    log("Please load dark, flat, and dark for flat SPE files",1)
    dark = np.zeros(img[0].shape)
    flat = np.ones(img[0].shape)
    #Select stars:
    selectstars()
    #spe.close() #In real version we'll close spe
    win.setupApsizeMenu()
        

#Determine the exposuretime of a SPE file without a footer
def getexptime(thisspe):
    #Input open SPE file
    #don't read lots of frames in large files
    numtoread = min([thisspe.get_num_frames(),11])
    tstamps = np.zeros(numtoread)
    for f in range(numtoread): 
        tstamps[f] = spe.get_frame(f)[1]['time_stamp_exposure_started']
    timediff = tstamps[1:numtoread]-tstamps[:numtoread-1]
    return np.round(np.median(timediff/1e6))


#Define all the stuff that needs to be done to each incoming frame
def processframe(i=0):
    global img,displayimg,rawtimes,backmed,backvar,framenum
    (thisframe,thistime) = spe.get_frame(i)
    #calibrate (doesn't do anything if calibration frames are not available):
    if darkExists: thisframe=(thisframe-dark)
    if flatExists: thisframe=thisframe/flat
    #read in frame
    img=np.transpose(thisframe)
    backgroundmed,backgroundvar=charbackground()
    #append stuff to global variables
    #Replace if this frame already exists, otherwise append
    if i <= framenum: #replace
        #log('Re-processing frame '+str(i)+' of '+str(framenum))
        rawtimes[i]=thistime['time_stamp_exposure_started']
        backmed[i]=backgroundmed
        backvar[i]=backgroundvar
    else: #append
        #log('Processing frame '+str(i)+' of '+str(framenum))
        rawtimes.append(thistime['time_stamp_exposure_started'])
        backmed.append(backgroundmed)
        backvar.append(backgroundvar)
    
    #make display image
    newdisplayimg=np.copy(img)
    newdisplayimg[0,0]=0
    imgvals = newdisplayimg.flatten()
    img99percentile = np.percentile(imgvals,99)
    newdisplayimg[newdisplayimg > img99percentile] = img99percentile
    #log("Framenum: "+str(framenum),2)
    #Replace if this frame already exists, otherwise append
    displayimg=newdisplayimg
    framenum=i
    
#Function to characterize the background to find stellar centroids accurately
#This should be done for each frame as it's read in
def charbackground():
    """Characterize the image background, median and variance

    for frame currenly held in img  
    """
    backgroundmed = biweight_location(img)
    backgroundvar = biweight_midvariance(img)
    return backgroundmed, backgroundvar        
 
#show the image to the widget
def displayFrame(autoscale=False,markstars=True):
    """Display an RBG image
    
    i is index to display
    Autoscale optional.
    Return nothing.
    """
    #Make sure i is in range
    if autoscale:
        #lowlevel=np.min(thisimg[thisimg > 0])
        lowlevel=np.percentile(displayimg[displayimg > 0],3)
        highlevel=np.max(displayimg)-1
        w5.setImage(np.array(displayimg),autoRange=True,levels=[lowlevel,highlevel],)
    else:
        w5.setImage(np.array(displayimg),autoRange=False,autoLevels=False)
    #Draw position circles:
    if markstars and len(stars) > 0:
        targs.setData([p[0] for p in stars[framenum]],[p[1] for p in stars[framenum]])
        targs.setSize(2.*apsizes[apsizeindex])
        targs.setPen(pencolors[0:numstars])
        
        
def selectstars():
    '''Select stars in the current frame.
    
    Click to select any number in the first image.
    Click to select numstars in later images to get following back on track.
    '''
    global selectingstars
    selectingstars = True
    

def improvecoords(x,y,i=framenum,pixdist=pixdist,fwhm=8.0,sigma=5.):
    """Improve stellar centroid position from guess value. (one at a time)

    #return the adjustment than needs to be made in x and y directions
    #NOTE: This may have trouble near edges.
    """
    #x=(1024/binning)-x
    #y=(1024/binning)-y
    #Keep track of motion
    delta = np.zeros(2)
    #Get image subregion around guess position
    subdata=img[x-pixdist:x+pixdist,y-pixdist:y+pixdist]
    #print subdata.shape
    sources = daofind(subdata - backmed[i], sigma*backvar[i], fwhm,
                      sharplo=0.1, sharphi=1.5, roundlo=-2.0, roundhi=2.0)
    #From what I can tell, daofind returns x and y swapped, so fix it
    returnedx = sources['ycentroid']
    returnedy = sources['xcentroid']
    
    if len(sources) != 0:
        strongsignal= np.argmax(sources['peak'])
        delta[0]+=returnedx[strongsignal]-pixdist
        delta[1]+=returnedy[strongsignal]-pixdist


    #check that unique source found
    '''
    if len(sources) == 0:
        log("Frame #"+str(i),1)
        log("WARNING: no sources found in searched region near ({:.2f}, {:.2f}).".format(x,y))
        #delta = [0,0] in this case
    else:
        if len(sources) > 1:
            log("Frame #"+str(i),1)
            log("WARNING: non-unique solution found for target near ({:.2f}, {:.2f}).".format(x,y))
            log(str(len(sources))+" signals in window.  Using brightest.")
        #Take brightest star found
    '''

    #handle stars that were not found #Move this outside this function
    """
    if [0,0] in delta and follow:
        meandeltax=np.mean(delta[np.where(delta[:,0] != 0),0])
        meandeltay=np.mean(delta[np.where(delta[:,1] != 0),1])
        delta[np.where(delta[:,0] == 0)] += [meandeltax,meandeltay]
    """

    return delta[0],delta[1]
















#### STAGE 2 ####

#Aperture details (provide a way to change these!)
apsizes=np.arange(1,11)
apsizeindex=3
r_in = 16.  #inner sky annulus radius #change in terms of binning eventually
r_out = 24. #outer sky annulus radius #change in terms of binning eventually
def setApSize(size):
    global apsizeindex
    log("Aperture size set to "+str(size)+" pixels.",1)
    #log("(Updates on next frame.)")
    if size in apsizes:
        apsizeindex=np.where(apsizes == size)[0][0]
        targs.setSize(2*size)# Currently doesn't update until next click/frame
        if stage > 1: 
            updatelcs(i=framenum)

compstar = 1 #which star to divide by
def setCompStar(s):
    global compstar
    compstar = s
    log("Now dividing by comparsion star #"+str(s),1)
    updatelcs(framenum)


#Phot results: variables to hold light curves and uncertainties
photresults=np.array([])


#Run the stage 2 loop
def stage2():
    global stars, spe, stage, hasFooter
    stagechange(2)
    #Add plot items for raw counts panel to plot
    for splot in rawcounts: w7.addItem(splot)
    #Make stars array an array of arrays of star coord arrays (yikes)
    # i.e, it needs to get pushed a level deeper
    stars=[stars]
    #Run photometry on the first frame
    dophot(0)
    updatelcs(i=0)
    updatehack()
    #Start timer that looks for new data
    timer2.start(exptime*1000.)
    timer3.start(10.*60*1000)#update every 10 minutes
    #This currently freezes up the UI.  Need to thread, but not enough time
    #to implement this currently.  Use a hack for now
    '''
    #Run the loop:
    fsize_spe_old = 0
    while not hasFooter:
        #Update only if there's new data
        fsize_spe_new = os.path.getsize(spefile)
        if fsize_spe_new > fsize_spe_old:
            spe = read_spe.File(spefile)
            numframes = spe.get_num_frames()
            log('Processing frames '+str(framenum)+'-'+str(numframes),1)
            while framenum < numframes:
                nextframe()
            if hasattr(spe, 'footer_metadata'): 
                hasFooter = True
                log('SPE footer detected. Data acquisition complete.',2)
                stagechange(3)
            spe.close()
        fsize_spe_old = fsize_spe_new
    '''


def updatehack():
    global spe, hasFooter, numframes
    #Only look for new data if not currently processing new data
    if not timer.isActive():
        fsize_spe_old = 0
        #Update only if there's new data
        fsize_spe_new = os.path.getsize(spefile)
        
        if fsize_spe_new > fsize_spe_old and stage ==2:
            spe = read_spe.File(spefile)
            numframes = spe.get_num_frames()
            if framenum+1==numframes-1:log('Processing frame '+str(framenum+1))
            else: log('Processing frames '+str(framenum+1)+'-'+str(numframes),1)        
            timer.start(100)
            #Update plots
            updatelcs(i=framenum)
            if hasattr(spe, 'footer_metadata'): 
                hasFooter = True
                log('SPE footer detected. Data acquisition complete.',2)
                stagechange(3)
                displayFrame()
        fsize_spe_old = fsize_spe_new

def nextframehack():
    #call nextframe until you're caught up
    global framenum,spe
    nextframe()
    updatelcs(i=framenum)
    if framenum >= numframes-1:
        timer.stop()
        updateft(i=framenum)
        spe.close()
        if stage==3:
            log("Image processing complete",2)
            timer3.stop()
            enddisplay()

#This timer catches up on photometry
timer = pg.QtCore.QTimer()#set up timer to avoid while loop
timer.timeout.connect(nextframehack)

#This timer checks for new data
timer2 = pg.QtCore.QTimer()
timer2.timeout.connect(updatehack)



#For demo purposes, read in the next frame of the spe file each time this is called
def nextframe():
    global stars
    #if stage == 2:
    oldcoords = stars[framenum]
    processframe(i=framenum+1) #Frame num increases here.
    newcoords=[]
    for coord in oldcoords:
        dx,dy = improvecoords(coord[0],coord[1],i=framenum)
        newcoords.append([np.floor(coord[0])+.5+dx,np.floor(coord[1])+.5+dy])        
    stars.append(newcoords)
    #Show the frame
    displayFrame(markstars=True)
    #Perform photometry
    dophot(i=framenum)
    #Update light curves
    #updatelcs(i=framenum) #only after all the new photometry is done.
    


    

def dophot(i):
    '''Do photometric measurements.
    
    Stars have been selected.  Do aperture photometry on given frame
    '''
    global photresults
    #print "dophot(i) called with i="+str(i)
    #Do the aperture photometry

    #The aperture_photometry() function can do many stars at once
    #But you must first do a background subtraction


    #We're going to save a lot of information in this step:
    #Total counts and uncertainty for every aperture size for every star
    #And eventually for every frame...

    #Note that the photometry package seems to reference x and y coords
    #as the tranpose of what we've been using.  Switch the order here:
    coords = [star[::-1] for star in stars[i]]
    thisphotometry = np.zeros((len(coords),len(apsizes)))
    for n in range(numstars):
        #Loop through the stars in the image
        #annulus_aperture = CircularAnnulus(coords[n], r_in=r_in, r_out=r_out)
        #print aperture_photometry(img[i],annulus_aperture).keys()
        #background_mean = aperture_photometry(img[i],annulus_aperture)['aperture_sum'][0]/annulus_aperture.area() 
        #NOTE on the above line: This should really be a median!
        #Issue 161 on photutils https://github.com/astropy/photutils/issues/161 is open as of 09/28/15
        gain = 12.63 #? From PI Certificate of Performance for "traditional 5MHz gain."  Confirm this value!
        #loop through aperture sizes
        for j,size in enumerate(apsizes):
            aperture = CircularAperture(np.array(coords[n]), r=size) 
            #phot = aperture_photometry(x-background_mean,aperture,error=backgroundvar,gain=gain)
            #Why am I getting negative numbers?
            #phot = aperture_photometry(img[i]-np.median(img),aperture)
            phot = aperture_photometry(img-backmed[i],aperture)
            thisphotometry[n,j]=phot['aperture_sum'][0]
    #print "photometry ",thisphotometry
    if i == 0:
        photresults = np.array([thisphotometry])
    else:
        #print "photresults dimensions are "+str(photresults.shape)
        #print "trying to append shape "+str(thisphotometry.shape)
        photresults = np.append(photresults,[thisphotometry],axis=0)
    #print "photresults dimensions are "+str(photresults.shape)
    #yay.  This deserves to all be checked very carefully, especially since the gain only affects uncertainty and not overall counts.


#define smoothing parameters For smoothed light curve for scipy.signal.savgol_filter
#winsize = 11.#win size in frames (must be odd)
sigma=3.
    

def enddisplay():
    #Change the end state of the display to a prettier image for the poster screenshot.
    global spe
    spe= read_spe.File(spefile)
    processframe(i=100)
    displayFrame(markstars=True)

#Update display.
def updatelcs(i):
    #Identify which points to include/exclude, up to frame i
    goodmask=np.ones(i+1, np.bool)
    goodmask[bad] = False
    badmask = np.zeros(i+1, np.bool)
    badmask[bad] = True
    targdivided = photresults[:i+1,0,apsizeindex]/photresults[:i+1,compstar,apsizeindex]
    times = np.arange(i+1)#Multiply by exptime for timestamps
    goodfluxnorm=targdivided[goodmask[:i+1]]/np.abs(np.mean(targdivided[goodmask[:i+1]]))
    s1.setData(exptime*times[goodmask[:i+1]],goodfluxnorm)
    #s2.setData(times[badmask[:i]],targdivided[badmask[:i]])
    l1.setData(exptime*times[goodmask[:i+1]],goodfluxnorm)
    #sl1.setData(times[goodmask[:i]],fluxsmoothed[goodmask[:i]])
    #Raw Counts:
    for j,splot in enumerate(rawcounts): splot.setData(exptime*times,photresults[:,j,apsizeindex])
    #Sky brightness
    sky.setData(exptime*times,backmed)

def updateftfromtimer():
    updateft(i=framenum)

def updateft(i=framenum):
        oversample=10. #Oversampling factor
        goodmask=np.ones(i+1, np.bool)
        goodmask[bad] = False
        targdivided = photresults[:i+1,0,apsizeindex]/photresults[:i+1,compstar,apsizeindex]
        goodfluxnorm=targdivided[goodmask[:i+1]]/np.abs(np.mean(targdivided[goodmask[:i+1]]))
        times = np.arange(i+1)#Multiply by exptime for timestamps
        #Fourier Transform   and smoothed lc  
        if goodmask.sum() > 2:
            #This all requires at least two points
            #Only update once per file read-in
            interped = interp1d(exptime*times[goodmask[:i+1]],goodfluxnorm-1.)
            xnew = np.arange(exptime*min(times[goodmask[:i]]),exptime*max(times[goodmask[:i+1]]),exptime)
            ynew = interped(xnew)
            #calculate FT
            amp = 2.*np.abs(fft(ynew,n=len(ynew)*oversample))#FFT
            amp /= float(len(ynew))
            freq = fftfreq(len(amp),d=exptime)
            pos = freq>=0 # keep positive part
            
            ft.setData(1e6*freq[pos],1e3*amp[pos])
            #Smoothed LC
            fluxsmoothed=filters.gaussian_filter1d(ynew[::oversample],sigma=sigma)
            ss1.setData(xnew[::oversample],fluxsmoothed)
    


def updateft_old(i=framenum): #update ft and smoothed lc
    if framenum < numframes:
        oversample=4. #Oversampling factor
        goodmask=np.ones(i+1, np.bool)
        goodmask[bad] = False
        targdivided = photresults[:i+1,0,apsizeindex]/photresults[:i+1,compstar,apsizeindex]
        goodfluxnorm=targdivided[goodmask[:i+1]]/np.abs(np.mean(targdivided[goodmask[:i+1]]))
        times = np.arange(i+1)#Multiply by exptime for timestamps
        #Fourier Transform   and smoothed lc  
        if goodmask.sum() > 2:
            #This all requires at least two points
            #Only update once per file read-in
            interped = interp1d(exptime*times[goodmask[:i+1]],goodfluxnorm-1.)
            xnew = np.arange(exptime*min(times[goodmask[:i+1]]),exptime*max(times[goodmask[:i+1]]),exptime/oversample)
            ynew = interped(xnew)
            #number of samples must be even for FT_continuous
            f=0
            H=0        
            if xnew.size % 2 == 0:
                f,H = FT_continuous(xnew,ynew)
            else:
                f,H = FT_continuous(xnew[:-1],ynew[:-1])
            H=2*np.sqrt(H.real**2 + H.imag**2.)/len(ynew)
            ft.setData(1e6*f[len(f)/2.:],1e3*H[len(f)/2.:])
            #Smoothed LC
            fluxsmoothed=filters.gaussian_filter1d(ynew[::oversample],sigma=sigma)
            ss1.setData(xnew[::oversample],fluxsmoothed)
            #QtGui.QApplication.processEvents()
        


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



## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        if len(sys.argv) > 1:
            defaultdir = sys.argv[1]
        QtGui.QApplication.instance().exec_()
