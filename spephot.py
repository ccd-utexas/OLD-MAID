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
from scipy import stats
from astropy.io import fits
from astropy.stats import biweight_location, biweight_midvariance
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
    dark=[]
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
            warnings+="No XML footer metadata. "
        
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
            warnings+="Shutter not closed for dark frame. "
        hdulist.close()
    else: warnings+="Invalid file type (must be SPE or FITS). "
    
    return dark, darkExp, warnings
    
    

#Load Flat frames
def openFlat(fname,darkForFlat,darkForFlatExp):
    '''
    Reduces, saves and returns flat from filename and dark data.
    
    Returns: 
        flat (2D numpy array)
        warnings (string)
    '''
    fname = str(fname)
    
    #Variable for master dark data
    flat=[]
    #Also collect any warning messages if there are any.
    warnings = ''
        
    if fname[-4:]=='.spe':
        if darkForFlat == []:
            warnings+="Import dark for reducting flats before importing flat SPE file. "
        else:
            fspe = read_spe.File(fname)
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
            
            #Write out fits file
            #Set up header
            prihdr = fits.Header()
            prihdr['OBJECT'] = 'flat'
            prihdr['IMAGETYP'] = 'flat'
            
            if hasattr(fspe, 'footer_metadata'):
                footer_metadata = BeautifulSoup(fspe.footer_metadata, "xml")
                ts_begin = footer_metadata.find(name='TimeStamp', event='ExposureStarted').attrs['absoluteTime']
                dt_begin = dateutil.parser.parse(ts_begin)
                prihdr['TICKRATE'] = int(footer_metadata.find(name='TimeStamp', event='ExposureStarted').attrs['resolution'])
                prihdr['DATE-OBS'] = str(dt_begin.isoformat())
                prihdr['XBINNING'] = footer_metadata.find(name="SensorMapping").attrs['xBinning']
                prihdr['YBINNING'] = footer_metadata.find(name="SensorMapping").attrs['yBinning']
                prihdr['INSTRUME'] = footer_metadata.find(name="Camera").attrs['model']
                prihdr['TRIGGER'] = footer_metadata.find(name='TriggerResponse').text
                prihdr['MODE'] = 1 #normalized
                prihdr['COMMENT'] = "SPE file has footer metadata"
                prihdr['EXPTIME'] = str(float(footer_metadata.find(name='ExposureTime').text)/1000.)
                flatexptime = np.round(float(footer_metadata.find(name='ExposureTime').text)/1000.)
                #check that dark exp time matches flat
                if flatexptime != darkForFlatExp:
                    warnings+="Exp times for dark and flat do not match! "
                    if darkForFlatExp == 0:
                       warnings+="Bias being used for flat subtraction. "
                #prihdr['SOFTWARE'] = footer_metadata.find(name='Origin')
                prihdr['SHUTTER'] = footer_metadata.find(name='Mode').text
                prihdr['REDUCED'] = dt.datetime.now().isoformat()
            else:
                prihdr['WARNING'] = "No XML footer metadata."
                warnings+="No XML footer metadata. "
                        #Set up fits object
            #Only write flat if properly dark subtracted:
            if warnings == '':
                hdu = fits.PrimaryHDU(flat,header=prihdr)
                flatpath = os.path.dirname(fname)
                fitsfilename = 'master_'+os.path.basename(fname).split('.spe')[0]+'.fits'
                hdu.writeto(os.path.join(flatpath, fitsfilename),clobber=True)
            #Close SPE
            fspe.close()
    #Option to load as Fits
    elif fname[-5:]=='.fits':
        hdulist = fits.open(fname)
        prihdr = hdulist[0].header
        flat=hdulist[0].data
        flatmode=stats.mode(flat.flatten())[0][0]
        print flatmode
        if flatmode != 1:
            warnings+= "Flat not properly normalized. Mode: "+str(flatmode)

        hdulist.close()            

    else: warnings+="Invalid file type (must be SPE)."
    
    return flat,warnings

    
def charbackground(img):
    """Characterize the image background, median and variance
    """
    backgroundmed = biweight_location(img)
    backgroundvar = biweight_midvariance(img)
    return backgroundmed, backgroundvar      
    
    
#Define all the stuff that needs to be done to each incoming frame
def reduceframe(frame,dark,flat):
    """Reduce frame with given dark, flat"""
    if dark != []: frame=(frame-dark)
    if flat != []: frame=frame/flat
    #transpose for correct orientation
    img=np.transpose(frame)
    return img
    
    
    