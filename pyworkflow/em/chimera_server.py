#!/usr/bin/env xmipp_python



from __future__ import print_function
from multiprocessing.connection import Listener, Client
from VolumeData import Array_Grid_Data
from VolumeViewer import volume_from_grid_data
from numpy import array, identity, dot
from numpy.linalg import inv
import chimera
from time import sleep
from threading import Thread
from chimera import runCommand
#from time import gmtime, strftime
#from datetime import datetime
import sys
#import socket

class ChimeraServer:
    
    def __init__(self):
        #print 'init'
        address = ''
        port = int(sys.argv[1])

        self.authkey = 'test'
        self.listener = Listener((address, port), authkey=self.authkey)
        self.vol_conn = self.listener.accept()
        self.openVolume()
        # Add basic handlers
        chimera.triggers.addHandler(chimera.APPQUIT, self.onAppQuit, None)
        # Add sub-classes handlers
        self.addHandlers()

        self.initListenShowJ()

    def addHandlers(self):
        """ Override this methods to specify which Chimera triggers want
        to be handled.
        """
        chimera.triggers.addHandler(chimera.MOTION_STOP, self.onMotionStop, None)####

    def openVolume(self):
        try:
            while True:
                if self.vol_conn.poll():
                    
                    msg = self.vol_conn.recv()
                    self.answer(msg)
                    if msg == 'end':#if you dont break cicle volume is never shown
                        break

                else:
                    sleep(0.01)
        except EOFError:
            print ('Lost connection to client')
            #should close app??

    def answer(self, msg):
        #print msg
        if msg == 'open_volume':
            data = self.vol_conn.recv()#objects are serialized by default
            #print data
            grid = Array_Grid_Data(data)
            self.volume = volume_from_grid_data(grid)

            
            #om = chimera.openModels
            #mlist = om.list()
            #m = mlist[0]
            #centerVec = m.bbox()[1].center().toVector()
            #xform = chimera.Xform.xform(1, 0, 0, -centerVec[0],
            #             0, 1, 0, -centerVec[1],
            #             0, 0, 1, -centerVec[2], orthogonalize=True)
            #m.openState.globalXform(xform)

            #runCommand("volume #0 step 1")


        elif msg == 'voxel_size':
            self.voxelSize = self.vol_conn.recv()
            cmd = "volume #0 voxelSize %s"%self.voxelSize
            #print cmd
            runCommand(cmd)
            #print reference spheres. Usefull for debuging, remove
            runCommand('shape sphere radius %s center %s,%s,%s color %s '%(self.voxelSize, -64 * self.voxelSize ,   0                  , 0 , 'orange'))
            runCommand('shape sphere radius %s center %s,%s,%s color %s '%(self.voxelSize,   0                  , -64 * self.voxelSize , 0 , 'forest green'))
            runCommand('shape sphere radius %s center %s,%s,%s color %s '%(self.voxelSize,   0                  ,   0                  , -64 * self.voxelSize  , 'magenta'))
            runCommand('shape sphere radius %s center %s,%s,%s color %s '%(self.voxelSize,  64 * self.voxelSize ,   0                  , 0 , 'red'))
            runCommand('shape sphere radius %s center %s,%s,%s color %s '%(self.voxelSize,   0                  ,  64 * self.voxelSize , 0 , 'green'))
            runCommand('shape sphere radius %s center %s,%s,%s color %s '%(self.voxelSize,   0                  ,   0                  , 64 * self.voxelSize , 'blue'))
            runCommand("focus")
            #end debug

        elif msg == 'command_list':
            commandList = self.vol_conn.recv()
            for command in commandList:
                runCommand(command)

    def initListenShowJ(self):

        self.listenShowJThread = Thread(target=self.listenShowJ)
        self.listenShowJThread.daemon = True
        self.listenShowJThread.start()

    def listenShowJ(self):
        try:
            om = chimera.openModels
            mlist = om.list()
            #assume that volume is 0 may be dangerous
            #this should be in init but it does not work there
            m = mlist[0]
            centerVec = m.bbox()[1].center().toVector()
            xform = chimera.Xform.xform(1, 0, 0, -centerVec[0],
                                        0, 1, 0, -centerVec[1],
                                        0, 0, 1, -centerVec[2], orthogonalize=True)
            m.openState.globalXform(xform)
            #end of this

            self.motion = identity(3)

            while True:
                if self.vol_conn.poll():

                    msg = self.vol_conn.recv()
                    #this rotate should not be generic
                    if msg == 'rotate':

                        matrix1 = self.vol_conn.recv()
                        matrix = dot(matrix1, inv(self.motion))#undo last rotation and put new one
                        self.motion = matrix1
                        xform = chimera.Xform.xform(matrix[0][0], matrix[0][1], matrix[0][2], 0,
                                       matrix[1][0], matrix[1][1], matrix[1][2], 0,
                                       matrix[2][0], matrix[2][1], matrix[2][2], 0, orthogonalize=True)

                        m.openState.globalXform(xform)
                    else:
                        sleep(0.01)
        except EOFError:
            print ('Lost connection to client')
            #should close app??

    
    def onMotionStop(self, trigger, extra, userdata):
        rx, ry , rz, t = self.volume.openState.xform.getCoordFrame()
        self.motion = array([[rx[0], ry[0], rz[0]], [rx[1], ry[1], rz[1]], [rx[2], ry[2], rz[2]]])
        printCmd('sending motion')
        self.vol_conn.send('motion_stop')
        self.vol_conn.send(self.motion)#send serialized motion
        printCmd('sended motion')

    def onAppQuit(self, trigger, extra, userdata):
        self.vol_conn.send('exit_server')
        self.listener.close()
        
            
def printCmd(cmd):
    pass
#            timeformat = "%S.%f"
#        print datetime.now().strftime(timeformat)
#        print cmd


class ChimeraVirusServer(ChimeraServer):

    def addHandlers(self):
        """ Override this methods to specify which Chimera triggers want
        to be handled.
        """
        pass
        #chimera.triggers.addHandler(chimera.MOTION_STOP, self.onMotionStop, None)####

    def answer(self, msg):
        """execute a single command and return values"""
        ChimeraServer.answer(msg)
        if msg == 'hk_icosahedron_lattice':
            from IcosahedralCage import cages
            h,k,radius,sym = self.vol_conn.recv()
            varray, tarray, hex_edges = cages.hk_icosahedron_lattice(h,k,radius,sym)
            self.vol_conn.send('hk_icosahedron_lattice')
            self.vol_conn.send(varray)

if len(sys.argv)> 1:
   serverName = sys.argv[2]
else:
   serverName = 'ChimeraServer'

serverClass = globals().get(serverName, None)
serverClass()

