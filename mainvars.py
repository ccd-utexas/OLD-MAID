# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 03:37:52 2016

This module holds all the top-level variables that everything needs access to.

@author: keatonb
"""

#imports
import numpy as np


stage=0 #start at 0
defaultdir = 'D:/sync_to_White_Dwarf_Archive/'#where to search for SPE files

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
darkExp=0 #exp time should match spe exptime
darkDark=False #shutter closed?
darkForFlat = []
darkForFlatExists=False
darkForFlatExp=0
darkForFlatDark=False
#Flat data
flat = []
flatExists=False
flatReduced=False #proper dark subtracted?
#Flag whether full reductions are being done (*correct* darks and flat)

#Phot results: variables to hold light curves and uncertainties
photresults=np.array([])

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
#List of median background counts:
backmed=[]
#List of background variances
backvar=[]
#Seeing for each star,frame:
seeing=[]
#Binning
binning=4

apsizes=np.arange(1,11)
apsizeindex=3

compstar = 1 #which star to divide by

fsize_spe_old = 0#Keep track if new spe file is larger that old one