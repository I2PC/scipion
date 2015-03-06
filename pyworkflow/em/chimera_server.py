#!/usr/bin/env xmipp_python



from multiprocessing.connection import Listener, Client
from VolumeData import Array_Grid_Data
from VolumeViewer import volume_from_grid_data
from numpy import array
import chimera
from time import sleep
from threading import Thread
from chimera import runCommand
from time import gmtime, strftime
from datetime import datetime
import sys
import socket


class ChimeraServer:
    
    def __init__(self):
        #print 'init'
        address = ''
        port = int(sys.argv[1])

        self.authkey = 'test'
        self.listener = Listener((address, port), authkey=self.authkey)
        self.vol_conn = self.listener.accept()
        self.openVolume()
        chimera.triggers.addHandler(chimera.MOTION_STOP, self.onMotionStop, None)
        chimera.triggers.addHandler(chimera.APPQUIT, self.onAppQuit, None)
        self.initListenShowJ()

    def openVolume(self):
        try:

            while True:
                if self.vol_conn.poll():
                    
                    msg = self.vol_conn.recv()
                    print msg
                    if msg == 'open_volume':
                        data = self.vol_conn.recv()#objects are serialized by default
                        #print data
                        grid = Array_Grid_Data(data)
                        self.volume = volume_from_grid_data(grid)
                        #runCommand("volume #0 step 1")
                        
                    elif msg == 'voxel_size':
                        voxelSize = self.vol_conn.recv()
                        cmd = "volume #0 voxelSize %s"%voxelSize
                        runCommand(cmd)
                        runCommand("focus")
                    
                    elif msg == 'draw_angular_distribution':
                        angulardist = self.vol_conn.recv()
                        for command in angulardist:
                            runCommand(command)

                    
                    elif msg == 'end':#if you dont break cicle volume is never shown
                        break
                else:
                    sleep(0.01)
        except EOFError:
            print 'Lost connection to client'
            #should close app??

    def initListenShowJ(self):

        self.listenShowJThread = Thread(target=self.listenShowJ)
        self.listenShowJThread.daemon = True
        self.listenShowJThread.start()

    def listenShowJ(self):
        try:

            while True:
                if self.vol_conn.poll():

                    msg = self.vol_conn.recv()
                    print msg
                    if msg == 'rotate':
                        matrix = self.vol_conn.recv()
                        print matrix
                        om = chimera.openModels
                        mlist = om.list()
                        m = mlist[0]
                        print type(matrix)
                        xform = chimera.Xform.xform(matrix[0][0],matrix[0][1], matrix[0][2], 0,
                                      matrix[1][0], matrix[1][1], matrix[1][2], 0,
                                      matrix[2][0], matrix[2][1], matrix[2][2], 0, orthogonalize=False)


                        m.openState.globalXform(xform)

                    else:
                        sleep(0.01)
        except EOFError:
            print 'Lost connection to client'
            #should close app??

    
    def onMotionStop(self, trigger, extra, userdata):
        print "Motion stop"
        rx, ry , rz, t = self.volume.openState.xform.getCoordFrame()
        self.motion = array([[rx[0], ry[0], rz[0]], [rx[1], ry[1], rz[1]], [rx[2], ry[2], rz[2]]])
        printCmd('sending motion')
        self.vol_conn.send('motion_stop')
        self.vol_conn.send(self.motion)#send serialized motion
        printCmd('sended motion')
            
            
            
            
    def onAppQuit(self, trigger, extra, userdata):
        print 'sended server exit'
        self.vol_conn.send('exit_server')
        self.listener.close()
        
            
def printCmd(cmd):
            timeformat = "%S.%f" 
#        print datetime.now().strftime(timeformat)
#        print cmd

ChimeraServer()
