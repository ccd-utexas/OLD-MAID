# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 05:50:43 2016

Define a custom process log with its own log method

@author: keatonb
"""

from pyqtgraph.Qt import QtGui


class ProcessLog(QtGui.QTextEdit):
    def __init__(self):
        super(ProcessLog,self).__init__()
        self.setReadOnly(True)
        
    def log(self,text,level=0):
        '''log messages to the process log and log file
        
        text is the message for the log
        level indicated how important it is:
        level=0: Routine background process: gray text;
        level=1: User influenced action: black text;
        level=2: Major change: bold black;
        level=3: Warning message: bold red; 
        '''
        text=str(text)
        colors = ['darkgray','black','black','red']
        prefix = ['','','','WARNING: ']
        fontweight = [50,50,75,75]
        if level in range(4):
            self.setTextColor(QtGui.QColor(colors[level]))
            self.setFontWeight(fontweight[level])
            self.append(prefix[level]+text)
        else: log('Level assigned to message "'+text+'" out of range.',level=3)
