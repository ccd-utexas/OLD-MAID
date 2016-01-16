# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 00:31:26 2016

Define a custom widget that calculates and plots the Fourier transform.

@author: keatonb
"""

#imports
import pyqtgraph as pg
import numpy as np
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq


class FTPlot(pg.PlotWidget):
    
    def __init__(self):
        super(FTPlot,self).__init__()
        self.ft = self.plot(pen='y') 
        self.setTitle("Fourier Transform")
        self.setLabel('left', 'amplitude (mma)')
        self.setLabel('bottom', 'freq (muHz)')
        
        #FT settings
        self.oversample=10. #Oversampling factor
        
    def plotft(self,freq,amp):
        self.ft.setData(freq,amp)

    def calcft(self,time,flux,exptime=1):
        #Only calculate the FT if there are >2 points
        if len(time) > 2:
        #interpolate over bad points #Might not be the best way to do this?
            interped = interp1d(time,flux)
            xnew = np.arange(min(time),max(time),exptime)
            ynew = interped(xnew)
            #calculate FT
            amp = 2.*np.abs(fft(ynew,n=len(ynew)*self.oversample))#FFT
            amp /= float(len(ynew))
            freq = fftfreq(len(amp),d=exptime)
            pos = freq>=0 # keep positive part
            
            self.plotft(1e6*freq[pos],1e3*amp[pos])
        
