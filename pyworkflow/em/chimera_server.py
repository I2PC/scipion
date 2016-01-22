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
    
    def __init__(self,centerVolume=True):
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
        self.centerVolume = centerVolume
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
                    if msg == 'open_volume':
                        data = self.vol_conn.recv()#objects are serialized by default
                        grid = Array_Grid_Data(data)
                        self.volume = volume_from_grid_data(grid)
                        
                    elif msg == 'voxel_size':
                        self.voxelSize = self.vol_conn.recv()
                        cmd = "volume #0 voxelSize %s"%self.voxelSize
                        runCommand(cmd)


                    elif msg == 'command_list':
                        commandList = self.vol_conn.recv()
                        for command in commandList:
                            runCommand(command)

                    elif msg == 'end':#if you don't break cicle volume is never shown
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
            self.centerVolume()


        elif msg == 'voxel_size':
            self.voxelSize = self.vol_conn.recv()
            cmd = "volume #0 voxelSize %s"%self.voxelSize
            #print cmd
            runCommand(cmd)
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

    def centerVolume(self):
        om = chimera.openModels
        mlist = om.list()
        #assume that volume is 0 may be dangerous
        #this should be in init but it does not work there
        self.model = mlist[0]
        centerVec = self.model.bbox()[1].center().toVector()
        xform = chimera.Xform.xform(1, 0, 0, -centerVec[0],
                                    0, 1, 0, -centerVec[1],
                                    0, 0, 1, -centerVec[2], orthogonalize=True)
        self.model.openState.globalXform(xform)
        #end of this

    def listenShowJ(self):
        try:
            while True:
                if self.vol_conn.poll():

                    msg = self.vol_conn.recv()
                    #this rotate should not be generic
                    if msg == 'rotate':

                        matrix1 = self.vol_conn.recv()
                        #undo last rotation and put new one. #Traslation is not undone, if user moves volume wrong translation applied
                        matrix = dot(matrix1, inv(self.rotation))

                        xform = chimera.Xform.xform(matrix[0][0], matrix[0][1], matrix[0][2], 0,
                                       matrix[1][0], matrix[1][1], matrix[1][2], 0,
                                       matrix[2][0], matrix[2][1], matrix[2][2], 0, orthogonalize=True)

                        self.model.openState.globalXform(xform)
                    else:
                        sleep(0.01)
        except EOFError:
            print ('Lost connection to client')
            #should close app??

    
    def onMotionStop(self, trigger, extra, userdata):
        rx, ry, rz, t = self.volume.openState.xform.getCoordFrame()
        self.rotation = array([[rx[0], ry[0], rz[0]], [rx[1], ry[1], rz[1]], [rx[2], ry[2], rz[2]]])

        self.vol_conn.send('motion_stop')
        self.vol_conn.send(self.rotation)#send serialized motion


    def onAppQuit(self, trigger, extra, userdata):
        self.vol_conn.send('exit_server')
        self.listener.close()
        
            
def printCmd(cmd):
    pass
#            timeformat = "%S.%f"
#        print datetime.now().strftime(timeformat)
#        print cmd


class ChimeraVirusServer(ChimeraServer):

    def __init__(self):
        ChimeraServer.__init__(self, centerVolume=True)

    def addHandlers(self):
        """ Override this methods to specify which Chimera triggers want
        to be handled.
        """
        pass
        #chimera.triggers.addHandler(chimera.MOTION_STOP, self.onMotionStop, None)####
        chimera.triggers.addHandler('selection changed', self.onselectionChanged, None)####

    def onselectionChanged(self, trigName, myData, trigData):
        sel = chimera.selection.currentGraphs()
        self.vol_conn.send('id')
        self.vol_conn.send(sel[0].id)#send serialized motion

    #not sure if  next function is answer or hanfleInit def handleInitMessage(self, msg):
    def answer(self, msg):
        """execute a single command and return values"""
        ChimeraServer.answer(msg)
        if msg == 'hk_icosahedron_lattice':
            from IcosahedralCage import cages
            h,k,radius,shellRadius,spheRadius,sym,sphere,color = \
                                      self.vol_conn.recv()
            ###
            #get vertexes of canonical triangle (20 per icosahedron)
            #get triangles defined by h k for canonical triangle
            corners, triangles, t_hex_edges = cages.hk_triangle(h, k)

            #get vertex for icosahedron
            #get vertex for each face
            from Icosahedron import icosahedron_geometry
            ivarray, itarray = icosahedron_geometry(sym)

            tlist = []
            #for a single face
            #map triangles to a single face in the given orientation
            for i0,i1,i2 in itarray:
                face = ivarray[i0], ivarray[i1], ivarray[i2]
                tmap = cages.triangle_map(corners, face)
                tlist.extend(cages.map_triangles(tmap, triangles))
                break#!!!!!!!!!!!!!!!!!!!!!!!!!
            va, ta = cages.surface_geometry(tlist, tolerance = 1e-5)

            from numpy import multiply
            multiply(va, shellRadius, va)    # Scale to requested radius            for point in va:
            for point in va:
                command = 'shape sphere radius %s center %s,%s,%s color %s '%\
                                                (spheRadius,
                                                 point[0],
                                                 point[1],
                                                 point[2],
                                                 color
                                                 )
                runCommand(command)
            self.vol_conn.send('axis')
            self.vol_conn.send(va)#send serialized motion
            ###!!!!hex_edges = array(t_hex_edges * len(itarray), intc)

            #show va spheres

if len(sys.argv)> 1:
   serverName = sys.argv[2]
else:
   serverName = 'ChimeraServer'

serverClass = globals().get(serverName, None)
serverClass()

