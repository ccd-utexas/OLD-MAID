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
        self.widthEdit.setMinimum(2)
        self.widthEdit.setMaximum(200)
        #self.widthEdit.setValue(currentwidth)
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
        dialog = KernelDialog(parent)
        result = dialog.exec_()
        kerneltype,width = dialog.kernelFormat()
        return (kerneltype,width, result == QtGui.QDialog.Accepted)


#set up a class that holds all the smoothing kernel information
class smoothingkernel:
    """Holds all smoothing kernel info"""
    kerneltype = 0
    width = 11 #points
    kernel=[]
    types = ['Uniform','Epanechnikov']
    def setkernel(self,kerneltype,width):
        if kerneltype == 1: #Epanechnikov
            u=(2.*np.arange(width)/(float(width)-1.))-0.5
            self.kernel = 0.75*(1.-u**2.)
            self.kernel /= np.sum(self.kernel)
            #processLog.log("Using "+self.types[kerneltype]+" smoothing kernel of width "+str(width))
        elif kerneltype == 0: #Uniform
            self.kernel = np.ones(width)/float(width)
            #processLog.log("Using "+self.types[kerneltype]+" smoothing kernel of width "+str(width))
    def openKernelDialog(self):
            dkerneltype,dwidth,daccepted = KernelDialog.getKernelFormat()
            if daccepted and (dkerneltype in range(len(types))) and (dwidth > 1): 
                self.setkernel(dkerneltype,dwidth)
    def __init__(self):
        self.setkernel(0,11)


class SmoothedLC(pg.PlotWidget):
    
    def __init__(self):
        super(SmoothedLC,self).__init__()
        #self.lcline = pg.PlotCurveItem()
        self.lcscatter = self.plot(brush=(255,0,0), pen='w',symbol='o') #Scatter plot
        self.setTitle("Smoothed Light Curve")
        self.setLabel('left', 'smoothed flux')
        self.setLabel('bottom', 'time (s)')
        
        self.kernel = smoothingkernel()


    def plotdata(self,time,flux,exptime):
        '''
        Smooth and plot light curves.
        '''
        #only update if there are enough points:
        if len(flux) > self.kernel.width:        
            #Fill bad points with nans:
            alltimes=np.arange(min(time),max(time)+exptime,exptime)/exptime
            filledflux = np.empty(alltimes.shape)
            filledflux[:] = np.NAN
            filledflux[(time/exptime).astype(int)] = flux
            #print type(filledflux)
            
            #convolve with smoothing kernel
            fluxsmoothed=convolve(filledflux,self.kernel.kernel,boundary='extend')
            self.lcscatter.setData(alltimes,fluxsmoothed)