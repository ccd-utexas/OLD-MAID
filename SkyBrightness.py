# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 19:45:12 2016

Plot the sky brightness over time.

@author: keatonb
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 13:20:23 2016

Define a custom widget that shows the divided light curve from raw data.

@author: keatonb
"""

import pyqtgraph as pg
from PyQt4 import QtCore
import numpy as np

class SkyBrightness(pg.PlotWidget):
    
    def __init__(self):
        super(SkyBrightness,self).__init__()
        #Vertical line shows currently displayed point
        self.currentind = pg.InfiniteLine(pen='y')
        self.addItem(self.currentind)
        #Scatter plot item
        self.skyscatter = self.plot(symbolBrush='w', symbolPen='w',pen='w',symbol='o', symbolSize=3) #Scatter plot
        self.setTitle("Sky Brightness")
        self.setLabel('left', 'median sky counts')
        self.setLabel('bottom', 'time (s)')
        
        self.numpts = 0
        
        #right click action
        #disable the normal context menu
        self.setMenuEnabled(False)
        self.scene().sigMouseClicked.connect(self.rightclicked)

    def plotsky(self,sky,exptime):
        '''
        Set the data
        
        - Sky is median counts per frame
        - exptime is exposure time
        '''
        self.sky=np.array(sky)
        self.numpts = self.sky.shape[0]
        self.exptime=exptime
        time=np.arange(self.numpts)*self.exptime

        #Plot result
        self.skyscatter.setData(time,self.sky)
        
        #Update vertical line marker if currently located at end of lc
        if self.currentind.value() >= max(time) - self.exptime:
            self.currentind.setValue(max(time))


    #Inspect location where right clicked
    def rightclicked(self,event):
        if event.button() == 2: #right click
            event.accept()
            x=self.skyscatter.getViewBox().mapSceneToView(event.scenePos()).x()
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


    def log(self,text,level=0):
        self.emit(QtCore.SIGNAL("log"),text,level)