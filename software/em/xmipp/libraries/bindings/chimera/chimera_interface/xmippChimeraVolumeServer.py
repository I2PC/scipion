# Server code ------------------------------
from multiprocessing.connection import Listener
import time
import numpy
import chimera
import SimpleSession
import pickle

class chimeraInterface():
    
    def __init__(self):
        '''Initialize variables, create socket, register events and open volume'''
        self.vProj = None
        self.address = ''
        self.port = 6000
        self.authkey = 'test'
        
        self.createSocket()
        self.registerEvent()
        self.openVolume()
    
    def createSocket(self):
        '''Create socket and wait for a new client'''
        print 'Waiting for client'
        self.listener = Listener((self.address, self.port), authkey=self.authkey) 
        self.remote_conn = self.listener.accept()
        print 'Got client ' + self.listener.last_accepted[0] + ':%d' % (self.listener.last_accepted[1])
          
         
    def openVolume(self):
        '''Wait for volume data and open in volume viewer'''
        try:
          while True:
            if self.remote_conn.poll():
              msg = self.remote_conn.recv()
              if type(msg) is numpy.ndarray:
                  from VolumeData import Array_Grid_Data
                  grid = Array_Grid_Data(msg)
              
                  from VolumeViewer import volume_from_grid_data
                  self.v = volume_from_grid_data(grid)
                  break
        
              #else:
                  #print 'msg: ' + msg
        
            else:
              time.sleep(0.01)
        except EOFError:
          print 'Lost connection to client'
          self.listener.close()
            
    #def motionStart(self, trigger, x, file):
        #print "Motion Start"
        
    def motionStop(self, trigger, x, file):
        '''Event executed when user unclick the volume''' 
        #print "Motion Stop"        
        rx, ry , rz, t = self.v.openState.xform.getCoordFrame()
        a = numpy.array([[rx[0], ry[0], rz[0]],
                 [rx[1], ry[1], rz[1]],
                 [rx[2], ry[2], rz[2]]])
        self.remote_conn.send(pickle.dumps(a))
                         
            
            
    def registerEvent(self):
        '''Register Events'''
        #chimera.triggers.addHandler(chimera.MOTION_START, self.motionStart, None)
        chimera.triggers.addHandler(chimera.MOTION_STOP, self.motionStop, None)

chimeraInterface()             
