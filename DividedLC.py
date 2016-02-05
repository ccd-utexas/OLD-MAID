# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 13:20:23 2016

Define a custom widget that shows the divided light curve from raw data.

@author: keatonb
"""

import pyqtgraph as pg
from PyQt4 import QtCore
import numpy as np

class DividedLC(pg.PlotWidget):
    
    def __init__(self):
        super(DividedLC,self).__init__()
        #Vertical line shows currently displayed point
        self.currentind = pg.InfiniteLine(pen='y')
        self.addItem(self.currentind)
        #Scatter plot item
        self.lcscatter = self.plot(symbolBrush='r', symbolPen='w',pen='w',symbol='o') #Scatter plot
        self.setTitle("Divided Light Curve")
        self.setLabel('left', 'rel. flux')
        self.setLabel('bottom', 'time (s)')
        
        #Keep track of bad points *not* to be plotted.
        self.badpoints=[]
        
        #Point Click action
        self.lcscatter.sigPointsClicked.connect(self.pointclicked)
        
        #right click action
        #disable the normal context menu
        self.setMenuEnabled(False)
        self.scene().sigMouseClicked.connect(self.rightclicked)

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
        self.numpts = self.flux.shape[0]
        self.time=np.arange(self.numpts)*self.exptime
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
        self.dividedlc=targdivided[goodmask]/np.abs(np.mean(targdivided[goodmask]))
        self.goodtime = self.time[goodmask]
        
        #Plot result
        self.lcscatter.setData(self.goodtime,self.dividedlc)
        
        #Update vertical line marker if currently located at end of lc
        if self.currentind.value() >= max(self.goodtime) - self.exptime:
            self.currentind.setValue(max(self.goodtime))
        
    #Flag point as bad if clicked on
    def pointclicked(self,item,points):
        #could be multiple overlapping points, so handle them individually
        for point in points:
            self.badpoints.append(point.pos()[0]/self.exptime)
            self.plotlc()
        
    #Inspect location where right clicked
    def rightclicked(self,event):
        if event.button() == 2: #right click
            event.accept()
            x=self.lcscatter.getViewBox().mapSceneToView(event.scenePos()).x()
            ind = np.round(x/self.exptime)
            if ind > self.numpts: #Follow newest frames
                ind = self.numpts -1
            elif ind < 0: ind = 0
            self.currentind.setValue(ind*self.exptime)
            #Emit a signal so other views can sync up with this            
            self.emit(QtCore.SIGNAL("inspectind"),ind)
            
    def setind(self,ind):
        if ind > self.numpts: #Follow newest frames
                ind = self.numpts -1
        elif ind < 0: ind = 0
        self.currentind.setValue(ind*self.exptime)
        
    def clearbadpoints(self):
        self.log("Deselecting "+str(len(self.badpoints))+" points.")
        self.badpoints = []
        
    def log(self,text,level=0):
        self.emit(QtCore.SIGNAL("log"),text,level)