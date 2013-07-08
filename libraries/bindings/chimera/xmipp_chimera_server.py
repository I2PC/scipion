#!/usr/bin/env xmipp_python


from pickle import dumps, loads
from multiprocessing.connection import Listener, Client
from VolumeData import Array_Grid_Data
from VolumeViewer import volume_from_grid_data
from numpy import array, ndarray
import chimera
from time import sleep
from threading import Thread
from os import system,environ
from chimera import runCommand
from time import gmtime, strftime
from datetime import datetime
import sys

class ChimeraServer:
    
    def __init__(self):
        #print 'init'
        self.address = ''
        #self.port = int(environ['XMIPP_CHIMERA_PORT'])
        self.port = 6000
        print sys.argv
        if len(sys.argv) > 2:#java case only
            self.port = int(sys.argv[2]) 
        elif len(sys.argv) > 1 and not 'xmipp_chimera_server.py' in sys.argv[1]:#avoid java case without port
            self.port = int(sys.argv[1]) 
        print self.port
        self.authkey = 'test'
        self.listener = Listener((self.address, self.port), authkey=self.authkey) 
        self.remote_conn = self.listener.accept()
        self.openVolume()
        
        chimera.triggers.addHandler(chimera.MOTION_STOP, self.onMotionStop, None)
        chimera.triggers.addHandler(chimera.APPQUIT, self.onAppQuit, None)
        self.initListen()
            
            
    def initListen(self):
        self.listen_thread = Thread(target=self.listen)
        self.listen_thread.daemon = True
        self.listen_thread.start()
    

    def listen(self):
        try:
            while True:
                if self.remote_conn.poll():
                    
                    msg = self.remote_conn.recv()
                    print msg
                    
                    if msg == 'exit_client':
                        break
                    else:
                        sleep(0.01)
        except EOFError:
            print 'Lost connection to client'
        finally:
            self.listener.close()
            runCommand('close all')
            runCommand('stop really')
            
            
            
    def openVolume(self):
        try:
            while True:
                if self.remote_conn.poll():
                    
                    msg = self.remote_conn.recv()
                    #print msg
                    if msg == 'open_volume':
                        data = self.remote_conn.recv()
                        #print 'volume data'
                        grid = Array_Grid_Data(data)
                        self.volume = volume_from_grid_data(grid)
                    if msg == 'draw_angular_distribution':
                        angulardist = self.remote_conn.recv()
                        for command in angulardist:
                            runCommand(command)
                    
                    if msg == 'end':    
                        break
                    else:
                        sleep(0.01)
        except EOFError:
            print 'Lost connection to client'
        

    
    def onMotionStop(self, trigger, extra, userdata):
        #print "Motion stop"
        rx, ry , rz, t = self.volume.openState.xform.getCoordFrame()
        self.motion = array([[rx[0], ry[0], rz[0]], [rx[1], ry[1], rz[1]], [rx[2], ry[2], rz[2]]])
        printCmd('sending motion')
        self.remote_conn.send('motion_stop')
        self.remote_conn.send(dumps(self.motion))#send serialized motion
        printCmd('sended motion')
            
            
            
            
    def onAppQuit(self, trigger, extra, userdata):

        self.remote_conn.send('exit_server')
        self.listener.close()
        
            
def printCmd(cmd):
            timeformat = "%S.%f" 
#        print datetime.now().strftime(timeformat)
#        print cmd

ChimeraServer()
