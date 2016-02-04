# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 13:20:23 2016

Define a custom widget that shows the divided light curve from raw data.

@author: keatonb
"""

import pyqtgraph as pg
import numpy as np

class DividedLC(pg.PlotWidget):
    
    def __init__(self):
        super(DividedLC,self).__init__()
        self.lcline = pg.PlotCurveItem()
        self.lcscatter = self.plot(brush=(255,0,0), pen='w',symbol='o') #Scatter plot
        self.setTitle("Divided Light Curve")
        self.setLabel('left', 'amplitude (mma)')
        self.setLabel('bottom', 'time (s)')
        
        #Keep track of bad points *not* to be plotted.
        self.badpoints=[]
        
        #Click action
        self.lcscatter.sigPointsClicked.connect(self.clicked)

    def setdata(self,flux,comp,apsizeindex,exptime):
        '''
        Set the data
        
        - Flux is a len(time)*N array for N stars (target is first).
        - Comp is the index of the comparison star.
        - apsize index is the index of the desired aperture.
        - exptime is exposure time in seconds.
        '''
        self.flux=flux
        self.comp=comp
        self.apsizeindex=apsizeindex
        self.exptime=exptime
        self.time=np.arange(self.flux.shape[0])*self.exptime
        self.plotlc()
    
    def plotlc(self):
        '''
        Plot the divided light curve.
        
        Don't plot points flagged as "bad"
        '''
        #Calculate divided light curve
        numpoints = len(self.time)
        goodmask=np.ones(numpoints, np.bool)
        goodmask[self.badpoints] = False
        compflux = self.flux[:,self.comp,self.apsizeindex]
        targdivided = self.flux[:,0,self.apsizeindex]/compflux
        dividedlc=targdivided[goodmask]/np.abs(np.mean(targdivided[goodmask]))
        time = self.time[goodmask]
        
        #Plot result
        self.lcline.setData(time,dividedlc)
        self.lcscatter.setData(time,dividedlc)
        
        
    def clicked(self,item,points):
        #could be multiple overlapping points, so handle them individually
        for point in points:
            self.badpoints.append(point.pos()[0]/self.exptime)
            self.plotlc()
        
    def clearbadpoints(self):
        self.badpoints = []