# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 03:48:58 2016

Define a custom widget that shows the smoothed light curve from the divided lc.

@author: keatonb
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
from astropy.convolution import convolve



#set up a dialog to change the kernel details:
class KernelDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(KernelDialog, self).__init__(parent)
        kerneltypes = ['Uniform','Epanechnikov']
        typeLabel = QtGui.QLabel("Kernel &type")
        self.typeEdit = QtGui.QComboBox()
        self.typeEdit.addItems(kerneltypes)
        #self.typeEdit.setCurrentIndex(currentind)
        typeLabel.setBuddy(self.typeEdit)
        widthLabel = QtGui.QLabel("Kernel &width")
        self.widthEdit = QtGui.QSpinBox()
        self.widthEdit.setMinimum(3)
        self.widthEdit.setMaximum(201)
        #self.widthEdit.setValue(width)
        widthLabel.setBuddy(self.widthEdit)
        self.buttons = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
                                        QtCore.Qt.Horizontal, self)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        
        grid = QtGui.QGridLayout(self)
        grid.addWidget(typeLabel,0,0)
        grid.addWidget(self.typeEdit,0,1)
        grid.addWidget(widthLabel,1,0)
        grid.addWidget(self.widthEdit,1,1)
        grid.addWidget(self.buttons, 3, 0)
        self.setLayout(grid)

        self.setWindowTitle("Define Smoothing Kernel")
    def kernelFormat(self):
        kerneltype=int(self.typeEdit.currentIndex())
        width=int(self.widthEdit.value())
        return (kerneltype,width)
        
    @staticmethod
    def getKernelFormat(parent = None):
        dialog = KernelDialog(parent,width=width)
        result = dialog.exec_()
        kerneltype,width = dialog.kernelFormat()
        return (kerneltype,width, result == QtGui.QDialog.Accepted)


#set up a class that holds all the smoothing kernel information
class smoothingkernel:
    """Holds all smoothing kernel info"""
    def __init__(self):
        self.kerneltype = 0
        self.width = 11 #points
        self.kernel=[]
        self.types = ['Uniform','Epanechnikov']
        self.setkernel(0,11)
    def setkernel(self,kerneltype,width):
        if kerneltype == 1: #Epanechnikov
            u=(2.*np.arange(width)/(float(width)-1.))-0.5
            self.kernel = 0.75*(1.-u**2.)
            self.kernel /= np.sum(self.kernel)
            self.log("Using "+self.types[kerneltype]+" smoothing kernel of width "+str(width))
        elif kerneltype == 0: #Uniform
            self.kernel = np.ones(width)/float(width)
            #processLog.log("Using "+self.types[kerneltype]+" smoothing kernel of width "+str(width))
    def openKernelDialog(self):
            dkerneltype,dwidth,daccepted = KernelDialog.getKernelFormat()
            if daccepted and (dkerneltype in range(len(self.types))) and (dwidth > 1): 
                self.setkernel(dkerneltype,dwidth)
    def log(self,text,level=0):
        self.emit(QtCore.SIGNAL("log"),text,level)


class SmoothedLC(pg.PlotWidget):
    
    def __init__(self):
        super(SmoothedLC,self).__init__()
        #Vertical line shows currently displayed point
        self.currentind = pg.InfiniteLine(pen='y')
        self.addItem(self.currentind)
        #Scatter plot shows light curve
        self.lcscatter = self.plot(brush=(255,0,0), pen='w',symbol='o') #Scatter plot
        self.setTitle("Smoothed Light Curve")
        self.setLabel('left', 'smoothed flux')
        self.setLabel('bottom', 'time (s)')
        
        self.numpts=0
        
        self.kernel = smoothingkernel()
        #self.connect(self.kernel,QtCore.SIGNAL("log"),self.log)
                       
        #right click action
        #disable the normal context menu
        self.setMenuEnabled(False)
        self.scene().sigMouseClicked.connect(self.rightclicked)


    def plotdata(self,time,flux,exptime):
        '''
        Smooth and plot light curves.
        '''
        self.exptime = exptime
        #only update if there are enough points:
        if len(flux) > self.kernel.width: 
            #Fill bad points with nans:
            alltimes=np.arange(min(time),max(time)+exptime,exptime)
            filledflux = np.empty(alltimes.shape)
            filledflux[:] = np.NAN
            filledflux[(((time-min(time)))/exptime).astype(int)] = flux
            #print type(filledflux)
            
            #convolve with smoothing kernel
            fluxsmoothed=convolve(filledflux,self.kernel.kernel,boundary='extend')
            self.lcscatter.setData(alltimes,fluxsmoothed)
            
            #Update vertical line marker if currently located at end of lc
            if self.currentind.value() >= self.numpts*self.exptime - self.exptime:
                self.currentind.setValue(max(alltimes))
            
            #update numpts
            self.numpts = max(alltimes)/exptime
            
    #Inspect location where right clicked
    def rightclicked(self,event):
        if event.button() == 2: #right click
            event.accept()
            x=self.lcscatter.getViewBox().mapSceneToView(event.scenePos()).x()
            ind = np.round(x/self.exptime)
            if ind > self.numpts: #Follow newest frames
                ind = self.numpts
            elif ind < 0: ind = 0
            self.currentind.setValue(ind*self.exptime)
            #Emit a signal so other views can sync up with this            
            #self.emit(QtCore.SIGNAL("inspectind"),ind)
        
    def log(self,text,level=0):
        self.emit(QtCore.SIGNAL("log"),text,level)