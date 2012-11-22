#!/usr/bin/env xmipp_python
'''
Created on Oct 22, 2012

@author: adrian
'''

import sys
import os
from subprocess import Popen, PIPE, STDOUT
import xmipp
import threading
from protlib_xmipp import getImageData
from numpy import ndarray, array, flipud
import pickle
import Queue


class openVolumesWithParameters:
    
    def __init__(self, inputFileName):
        '''Initialize variables and call openVolume'''
        self.projectionImage = None
        self.image = None
        self.file = inputFileName
        self.address = ''
        self.port = 6000
        self.authkey = 'test'
        #self.isChange = True
        self.openVolume()
           
    def openVolume(self):
        '''Open Volume, send through the socket, create first projection and visualize it and create parallel thread for incoming events'''
        self.image = xmipp.Image(self.file)
        xdim, ydim, zdim, n = self.image.getDimensions()
        Z = getImageData(self.image)
        self.image.convert2DataType(xmipp.DT_DOUBLE)
        
        from multiprocessing.connection import Client
        address = (self.address, self.port)
        self.conn = Client(address, authkey=self.authkey)
        self.conn.send(Z)
    
        self.projection = self.image.projectVolumeDouble(0,0,0)
            
        self.process_check = True
        self.process_thread = threading.Thread(target=self.communication)
        self.process_thread.start()
            
        from protlib_gui_figure import ImageWindow
        self.projectionImage = ImageWindow(image=self.projection,label="Projection")
        
        self.projectionImage.root.bind('<<NewProjection>>', self.update)
        #self.update()
        self.projectionImage.show()
            
        
    def communication(self):
        '''Function executed by parallel thread
        Wait for any data comming from the socket and calculate projection'''
        keep_running = True
        while keep_running:
            msg = self.conn.recv()
            
            msg2= pickle.loads(msg)
            b=array(msg2)
            
            rot1,tilt1,psi1 = xmipp.Euler_matrix2angles(b)
#            print "rot1,tilt1,psi1"
#            print rot1,tilt1,psi1
            
            self.projection = self.image.projectVolumeDouble(rot1,tilt1,psi1)
            #self.isChange = True
            self.projectionImage.root.event_generate('<<NewProjection>>', when='tail')
            
            
    def update(self, event):
        '''Update projection in viewer when an event is detected'''
        #print event
        #if self.isChange:
            #self.projection = self.image.projectVolumeDouble(rot, tilt, psi)
        Z = flipud(getImageData(self.projection))
        self.projectionImage.updateData(Z) 
        #self.isChange = False
        #self.projectionImage.root.after(50, self.update)                 
    
        

if __name__ == '__main__':
    def command_line_options():
        """ add command line options here"""
        import optparse
        _usage = "usage: %prog [options] Example:   %prog -i task.xml -q False"
        global operation
        parser = optparse.OptionParser(_usage)
        
        parser.add_option("-i", "--inputFileName", dest="inputFileName",
                          default="/dev/stdin", type="string",
                          help="File to display")                           
        
        (options, args) = parser.parse_args()
        return(options.inputFileName)
    
    inputFileName  = command_line_options()
    instance = openVolumesWithParameters(inputFileName)