#!/usr/bin/env xmipp_python


from pickle import dumps, loads
from multiprocessing.connection import Listener, Client
from VolumeData import Array_Grid_Data
from VolumeViewer import volume_from_grid_data
from numpy import array, ndarray
import chimera
from time import sleep

class ChimeraServer:
    
    def __init__(self):
        print 'init'
        self.volfile = file
        self.address = ''
        self.port = 6000
        self.authkey = 'test'
        self.listener = Listener((self.address, self.port), authkey=self.authkey) 
        self.remote_conn = self.listener.accept()
        self.quit = False
        self.listen()
        chimera.triggers.addHandler(chimera.MOTION_STOP, self.onMotionStop, None)
        chimera.triggers.addHandler(chimera.APPQUIT, self.onAppQuit, None)
    

    def listen(self):
        print 'waiting for volume data to open in volume viewer'
        try:
            while True:
                if self.remote_conn.poll():
                    msg = self.remote_conn.recv()
                    print 'connection received'
                    if type(msg) is ndarray:
                        grid = Array_Grid_Data(msg)
                        self.volume = volume_from_grid_data(grid)
                        break
                else:
                    sleep(0.01)
        except EOFError:
            print 'Lost connection to client'
            
    
    
    def onMotionStop(self, trigger, extra, userdata):
        print "Motion stop"
        rx, ry , rz, t = self.volume.openState.xform.getCoordFrame()
        self.motion = array([[rx[0], ry[0], rz[0]], [rx[1], ry[1], rz[1]], [rx[2], ry[2], rz[2]]])
        #print self.motion
        self.remote_conn.send(dumps(self.motion))#send serialized motion
            
    def onAppQuit(self, trigger, extra, userdata):
        self.quit = True
        self.remote_conn.send(dumps('exit'))
        self.listener.close()

ChimeraServer()