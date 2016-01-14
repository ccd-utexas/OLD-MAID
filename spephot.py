# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 02:59:13 2016

Collect all the functions related to frame reduction and photometry.

@author: keatonb
"""

#imports
import numpy as np
import os
import datetime as dt
import dateutil.parser
from astropy.io import fits
from bs4 import BeautifulSoup
import read_spe


def openDark(fname):
    '''
    Reduces, saves and returns dark from filename.
    
    Returns: 
        dark (2D numpy array)
        darkExp (int seconds)
        warnings (string)
    '''
    fname = str(fname)
    
    #Variable for master dark data
    dark=None    
    
    #Be ready to store the exposure time
    darkExp=np.nan
    
    #Also collect any warning messages if there are any.
    warnings = ''
        
    
    if fname[-4:]=='.spe':
        #processLog.log("Opening dark file "+fname,1)
        dspe = read_spe.File(fname)
        num_darks=dspe.get_num_frames()
        
        #get all frames in SPE file
        #stack as 3D numpy array
        (frames,_)=dspe.get_frame(0)
        frames=np.array([frames])
        for i in range(1,num_darks):
            (thisframe,_)=dspe.get_frame(i)
            frames=np.concatenate((frames,[thisframe]),0)
        #Calculate the master dark
        dark=np.median(frames,axis=0)
        
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
            darkExp=np.round(float(footer_metadata.find(name='ExposureTime').text)/1000.)
            prihdr['EXPTIME'] = str(float(footer_metadata.find(name='ExposureTime').text)/1000.)
            prihdr['SHUTTER'] = footer_metadata.find(name='Mode').text
            if footer_metadata.find(name='Mode').text != 'AlwaysClosed':
                prihdr['WARNING'] = 'Shutter not closed for dark frame.'
                warnings+="Shutter not closed for dark frame. "
        else:
            prihdr['WARNING'] = "No XML footer metadata."
            warnings+="No XML footer metadata."
        
        #Set up fits object
        hdu = fits.PrimaryHDU(dark,header=prihdr)
        darkpath = os.path.dirname(fname)
        fitsfilename = 'master_'+os.path.basename(fname).split('.spe')[0]+'.fits'
        hdu.writeto(os.path.join(darkpath, fitsfilename),clobber=True)
        #Close SPE
        dspe.close()
    
    #option to load as fits
    elif fname[-5:]=='.fits':
        hdulist = fits.open(fname)
        prihdr = hdulist[0].header
        dark=hdulist[0].data
        darkExp = np.round(float(prihdr['EXPTIME']))
        if prihdr['SHUTTER'] != 'AlwaysClosed':
            prihdr['WARNING'] = 'Shutter not closed for dark frame.'
            warnings+="Shutter not closed for dark frame."
        hdulist.close()
    else: warnings+="Invalid file type (must be SPE or FITS)."
    
    return dark, darkExp, warnings
    
    
    
    
    