# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 07:52:35 2015

@author: keatonb
"""

#!/usr/bin/python
# -*- coding: utf-8 -*-



# Imports.
# Standard libraries.
from __future__ import absolute_import, division, print_function
import os

#CD to script location
os.chdir('C:/Users/admin/Documents/GitHub/tsphot')

# More Standard libraries.
import sys
import time
import shutil
import tempfile
import traceback
import datetime as dt
# Installed packages.
from astropy.io import fits
from PyQt4 import QtGui, QtCore
# Local modules.
import read_spe


defaultdir = 'D:/sync_to_White_Dwarf_Archive/'#where to search for SPE files

class toAutoguider(QtGui.QWidget):
    
    def __init__(self):
        super(toAutoguider, self).__init__()
        self.initUI()
        
    def initUI(self):      

        self.log = QtGui.QTextEdit(self)
        self.log.move(20, 20)
        self.log.setReadOnly(True)
        
        self.btn = QtGui.QPushButton('New SPE File', self)
        self.btn.move(20, 220)
        self.btn.clicked.connect(self.showDialog)
        
        self.btn2 = QtGui.QPushButton('Close Program', self)
        self.btn2.move(145, 220)
        self.btn2.clicked.connect(QtGui.qApp.quit)
                
        
        self.setGeometry(300, 300, 290, 260)
        self.setWindowTitle('Sending Data to Autoguider')
        self.show()
        
    def showDialog(self):
        
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', 
                defaultdir)
        
        runloop(str(fname))
        #with f:        
        #    data = f.read()
        #    self.log.setText(data) 


# Make the App have a window and dock area.         
app = QtGui.QApplication([])
win = toAutoguider()
#win.setCentralWidget(area)
#win.resize(1500,800)
win.setWindowTitle('Send to Autoguider')

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
        win.log.setTextColor(QtGui.QColor(colors[level]))
        win.log.setFontWeight(fontweight[level])
        win.log.append(prefix[level]+text)
    else: log('Level assigned to message "'+text+'" out of range.',level=3)




# Define function to write last frame to FITS.
def write_last_frame_to_fits(fpath_spe, fpath_fits, pixscale=0.36, orient='UPPER_LEFT'):
    """Write the last frame of an SPE file to FITS format.

    Parameters
    ----------
    fpath_spe : string
        Absolute path to SPE file.
    fpath_fits : string
        Absolute path to temporary FITS file.
    pixscale : {0.36}, float, optional
        Plate scale (arcsec/superpixel).
        Default value 0.36 arcsec/superpix is for ProEM 1024B at Cassegrain focus on McDonald 2.1m.
    orient : {'UPPER_LEFT', 'LOWER_LEFT', 'UPPER_RIGHT', 'LOWER_RIGHT'}, string, optional
        Location of pixel (0, 0) to match frame orientation for North-up, East-left.
        Default value of (0, 0) in 'UPPER_LEFT' is for ProEM 2014B with default LightField settings.

    Returns
    -------
    None

    Notes
    -----
    - The default for DS9 is to display (0, 0) in lower left.
    - If DS9 has temp.fits open, DS9 locks the file from being written to.
    - Do not use the file creation time for your final reductions.
      Use the timestamp from the SPE metadata XML footer, which is created after clicking "Stop".
    """
    # Open SPE, load frame, check for existing FITS file.
    spe_ctime = os.path.getctime(fpath_spe)
    spe = read_spe.File(fpath_spe)
    #is there a footer:    
    nofooter = not hasattr(spe, 'footer_metadata')
    
    (frame, metadata) = spe.get_frame(-1)
    if os.path.isfile(fpath_fits):
        os.remove(fpath_fits)
    orient_options = ['UPPER_LEFT', 'UPPER_RIGHT', 'LOWER_LEFT', 'LOWER_RIGHT']
    if orient not in orient_options:
        raise IOError(("Invalid argument for orient:\n" +
                       "Given: orient = {orient}\n" +
                       "Valid options: {orient_options}").format(orient=orient, orient_options=orient_options))
    # Create FITS file.
    hdu = fits.PrimaryHDU(frame)
    hdulist = fits.HDUList([hdu])
    prihdr = hdulist[0].header
    fctime = dt.datetime.utcfromtimestamp(spe_ctime).isoformat()
    prihdr['FCTIME']   = (fctime,
                          'SPE file creation time (UTC, ISO 8601)')
    fctimeux = spe_ctime
    prihdr['FCTIMEUX'] = (fctimeux,
                          'SPE file creation time (Unix time)')
    framenum = metadata['frame_tracking_number']
    prihdr['FRAMENUM'] = (framenum,
                          'Frame tracking number. FRAMENUM >= 1')
    expstart = metadata['time_stamp_exposure_started'] / 1.0E6
    prihdr['EXPSTART'] = (expstart,
                          'Time (s) from "Run Inf"/"Acquire" to exp start')
    expend = metadata['time_stamp_exposure_ended'] / 1.0E6
    prihdr['EXPEND']   = (expend,
                          'Time (s) from "Run Inf"/"Acquire" to exp end')
    exptime = expend - expstart
    prihdr['EXPTIME']  = (exptime,
                          'Exposure time (s). EXPTIME = EXPEND - EXPSTART')
    prihdr['PIXSCALE'] = (pixscale,
                          'Plate scale (arcsec/superpixel)')
    prihdr['ORIENT']   = (orient,
                          'Location of pixel (0,0) to display N-up, E-left')
    # Write FITS file. Close FITS and SPE files.
    hdulist.writeto(fpath_fits)
    hdulist.close()
    spe.close()
    #return where there is a footer
    return nofooter
    


# Define sleep time and message.
sleep_time = .5 # seconds
sleep_msg = ("INFO: Sleeping for {num} seconds between iterations.").format(num=sleep_time)
log(sleep_msg)


# Continuously write the last SPE frame to a temporary FITS file.
# Note: If DS9 has temp.fits open, DS9 locks the file from being written to.
fpath_spe  = ""
fdir_fits = os.path.join(tempfile.gettempdir(), 'Guide82')
if not os.path.exists(fdir_fits):
    os.makedirs(fdir_fits)
fname_fits_writing = 'temp_WRITING.fits'
fpath_fits_writing = os.path.join(fdir_fits, fname_fits_writing)
fname_fits = 'temp.fits'
fpath_fits = os.path.join(fdir_fits, fname_fits)
fsize_spe_old = 0
write_fits = None

def runloop(fname):
    # Define file paths.
    global fpath_spe
    fpath_spe  = os.path.abspath(fname)
    if not os.path.isfile(fpath_spe):
        raise IOError(("File does not exist: {fpath_spe}").format(fpath_spe=fpath_spe))
    log(("SPE file:\n"+
           "{fs}\n"+
           "FITS file for writing:\n"+
           "{ffw}\n"+
           "FITS file:\n"+
           "{ff}").format(fs=fpath_spe, ffw=fpath_fits_writing, ff=fpath_fits),2)
        

    timer.start(1000.*sleep_time)

def loop():
    global fsize_spe_old, write_fits

    # Check the SPE file size for later comparison.
    fsize_spe_new = os.path.getsize(fpath_spe)
    # On the first iteration, always write out a new FITS file...
    if write_fits is None:
        write_fits = True
    # ...otherwise check that a new FITS file should be written out.
    else:
        # If LightField wrote to the SPE file, write out a new FITS file...
        if fsize_spe_old < fsize_spe_new:
            write_fits = True
        # ...otherwise do nothing.
        else:
            write_fits = False

    nofooter = None            
    # Write to a temporary file then move as an atomic operation to prevent
    # file access conflicts with Guide82.
    if write_fits:
        try:
            nofooter = write_last_frame_to_fits(fpath_spe=fpath_spe,
                                     fpath_fits=fpath_fits_writing)
            shutil.move(src=fpath_fits_writing, dst=fpath_fits)
            log(("INFO: FITS file updated at {dt_now}.\n"+
                   "FITS file:\n"+
                   "{fpath_fits}").format(dt_now=dt.datetime.now(),
                                            fpath_fits=fpath_fits))
        except (WindowsError, IOError) as err:
            log(err)
            log(("A program may have locked the FITS file.\n"+
                   "Close the program (e.g. DS9) to allow writing to the FITS file.\n"+
                   "FITS file was not updated:\n"+
                   "{fpath_fits}").format(fpath_fits=fpath_fits),3)
    # Set comparison variables for next iteration then sleep.
    fsize_spe_old = fsize_spe_new
    write_fits = False
    
    #end loop if there's a footer
    if nofooter == False:
        timer.stop()
        log("SPE footer detected.  Sending data to Autoguider stopped.",2)

timer = QtCore.QTimer()#set up timer to avoid while loop
timer.timeout.connect(loop)



## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        if len(sys.argv) > 1:
            runloop(sys.argv[1])#extra argument is spe filepath
        QtGui.QApplication.instance().exec_()