# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 21:02:30 2016

Define a custom image display with its own display methods

@author: keatonb
"""

#imports
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import numpy as np

class ImageDisplay(pg.ImageView):
    def __init__(self):
        super(ImageDisplay,self).__init__()
        self.ui.roiBtn.hide()
        #self.ui.normBtn.hide() #Causing trouble on Windows?
        
        #set up the details for creating display image:
        self.clippercent = 99 #Clip cosmic rays and such
        self.displaypercents = [5,100]
        
        #Setup star markers:
        self.setupMarkers()        
        
    def displayImage(self,img):
        displayimg=np.copy(img)
        #cut off high outliers: cosmic rays and such
        cutoff = np.percentile(displayimg,self.clippercent)
        displayimg[displayimg > cutoff] = cutoff
        #Set black/while levels
        levels = np.percentile(displayimg,self.displaypercents)
        #set one pixel to zero for scaling
        displayimg[0,0]=0
        self.setImage(displayimg,autoRange=False,levels=levels)
        
    def setupMarkers(self):
        '''Set up plot item for showing apertures'''
        #print QtGui.QColor.colorNames() for available names.
        stringcolors=['red','green','blue','magenta','orange','yellow',
                      'darkred','darkgreen','darkblue','darkmagenta','darkorange','darkgoldenrod',
                      'hotpink','seagreen','skyblue','salmon','brown','lightyellow']
        self.pencolors = [pg.mkPen(QtGui.QColor(c), width=3) for c in stringcolors]
        self.aps = pg.ScatterPlotItem(brush=None, pen=self.pencolors[0],symbol='o',pxMode=False,size=4)
        self.addItem(self.aps)
        
    def displayApertures(self,coords,apsize):
        '''draw the new apertures in position'''
        self.aps.setData([p[0] for p in coords],[p[1] for p in coords])
        self.aps.setSize(2.*apsize)
        self.aps.setPen(self.pencolors[0:len(coords)])