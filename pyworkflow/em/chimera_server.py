#!/usr/bin/env xmipp_python



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

        if len(sys.argv) == 3:
            clientport = int(sys.argv[2])
            self.serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.serversocket.bind((address, clientport))
            #become a server socket
            self.serversocket.listen(1)
            self.listen_thread = Thread(target=self.listenClients)
            self.listen_thread.daemon = True
            self.listen_thread.start()


    def listenClients(self):


        while True:
            print 'on listening'
            try:
                (clientsocket, address) = self.serversocket.accept()

                print "connection accepted"

                msg = clientsocket.recv(1024)#should be a single message, so no loop
                zrot = chimera.Xform.zRotation(30)	# rotation by 30 degrees
                om = chimera.openModels.list()[0]
                om.openState.globalXform(zrot)

                print 'msg %s' % msg
                clientsocket.sendall(msg)
                if msg == 'end':
                    break

                self.clientsocket.close()

            except EOFError:
                print 'Lost connection to client'




    def openVolume(self):
        try:

            while True:
                if self.vol_conn.poll():
                    
                    msg = self.vol_conn.recv()
                    #print msg
                    if msg == 'open_volume':
                        data = self.vol_conn.recv()#objects are serialized by default
                        #print data
                        print type(data)
                        grid = Array_Grid_Data(data)
                        self.volume = volume_from_grid_data(grid)
                        #runCommand("volume #0 step 1")
                        
                    elif msg == 'voxelSize':
                        voxelSize = self.vol_conn.recv()
                        cmd = "volume #0 voxelSize %s"%voxelSize
                        runCommand(cmd)
                        runCommand("focus")
                    
                    elif msg == 'draw_angular_distribution':
                        angulardist = self.vol_conn.recv()
                        for command in angulardist:
                            runCommand(command)
                    
                    elif msg == 'end':    
                        break
                else:
                    sleep(0.01)
        except EOFError:
            print 'Lost connection to client'
            #should close app??
        

    
    def onMotionStop(self, trigger, extra, userdata):
        #print "Motion stop"
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
