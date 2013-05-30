#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
 '''

from xmipp import *
from multiprocessing.connection import Client
from protlib_xmipp import getImageData
from protlib_gui_figure import ImageWindow
from pickle import dumps, loads
from threading import Thread
from os import system
from numpy import array, ndarray, flipud
from time import gmtime, strftime
from datetime import datetime
from os.path import exists
from decimal import *

class XmippChimeraClient:
    
    def __init__(self, volfile, angulardistfile=None, spheres_color='red', spheres_distance='default', spheres_maxradius='default'):
        
        if volfile is None:
            raise ValueError(volfile)
        if '@' in volfile:
            [index, file] = volfile.split('@'); 
        else :
            file = volfile
        if not exists(file):
            raise ValueError(file)
        if not angulardistfile is None:
            if not(exists(angulardistfile)):
                raise ValueError(angulardistfile)
        
            
        self.volfile = volfile
        self.angulardistfile = angulardistfile
        
        self.address = ''
        self.port = 6000
        self.authkey = 'test'
        self.client = Client((self.address, self.port), authkey=self.authkey)
        printCmd('initVolumeData')
        self.initVolumeData()
        self.spheres_color = spheres_color
        self.spheres_distance = float(spheres_distance) if not spheres_distance == 'default' else max(self.xdim, self.ydim, self.zdim)
#        print self.spheres_distance
        self.spheres_maxradius = float(spheres_maxradius) if not spheres_maxradius == 'default' else 0.02 * self.spheres_distance
      
        printCmd('openVolumeOnServer')
        self.openVolumeOnServer(self.vol)
    
    
    def loadAngularDist(self):
        md = MetaData(self.angulardistfile)
        maxweight = md.aggregateSingle(AGGR_MAX, MDL_WEIGHT)
        minweight = md.aggregateSingle(AGGR_MIN, MDL_WEIGHT)
        interval = maxweight - minweight
        if interval < 1:
            interval = 1
        minweight = minweight - 1#to avoid 0 on normalized weight
        self.angulardist = []  
        x2=self.xdim/2
        y2=self.ydim/2
        z2=self.zdim/2
        #cofr does not seem to work!
        #self.angulardist.append('cofr %d,%d,%d'%(x2,y2,z2))
        for id in md:
            
            rot = md.getValue(MDL_ANGLE_ROT, id)
            tilt = md.getValue(MDL_ANGLE_TILT, id)
            psi = md.getValue(MDL_ANGLE_PSI, id)
            weight = md.getValue(MDL_WEIGHT, id)
            weight = (weight - minweight)/interval

            x, y, z = Euler_direction(rot, tilt, psi)
            radius = weight * self.spheres_maxradius
            x = x * self.spheres_distance+x2
            y = y * self.spheres_distance+y2
            z = z * self.spheres_distance+z2
            command = 'shape sphere radius %s center %s,%s,%s color %s '%(radius, x, y, z, self.spheres_color)
            self.angulardist.append(command)    
            printCmd(command)
            
    def send(self, cmd, data):
#        print cmd
        self.client.send(cmd)
        self.client.send(data)
        
        
    def openVolumeOnServer(self, volume):
         self.send('open_volume', volume)
         if not self.angulardistfile is None:
             self.loadAngularDist()
             self.send('draw_angular_distribution', self.angulardist)
         self.client.send('end')
        

    def initListenThread(self):
            self.listen_thread = Thread(target=self.listen)
            self.listen_thread.daemon = True
            self.listen_thread.start()
         
    
    def listen(self):
        
        self.listen = True
        try:
            while self.listen:
                #print 'on client loop'
                msg = self.client.recv()
                self.answer(msg)
                    
        except EOFError:
            print 'Lost connection to server'
        finally:
            self.exit()
            
            
    def exit(self):
            self.client.close()#close connection


    def initVolumeData(self):
        self.image = Image(self.volfile)
        self.image.convert2DataType(DT_DOUBLE)
        self.xdim, self.ydim, self.zdim, self.n = self.image.getDimensions()
        printCmd("size %dx %dx %d"%(self.xdim, self.ydim, self.zdim))
        self.vol = getImageData(self.image)
        
    def answer(self, msg):
        if msg == 'exit_server':
            self.listen = False



class XmippProjectionExplorer(XmippChimeraClient):
    
    def __init__(self, volfile, angulardistfile=None, spheres_color='red', spheres_distance='default', spheres_maxradius='default', size='default', padding_factor=1, max_freq=0.5, spline_degree=2):

        XmippChimeraClient.__init__(self, volfile, angulardistfile, spheres_color, spheres_distance, spheres_maxradius)
        
        self.projection = Image()
        self.projection.setDataType(DT_DOUBLE)
        #0.5 ->  Niquiest frequency
        #2 -> bspline interpolation
        #print 'creating Fourier Projector'
        self.fourierprojector = FourierProjector(self.image, padding_factor, max_freq, spline_degree)
        self.fourierprojector.projectVolume(self.projection, 0, 0, 0)

        printCmd('initListenThread')
        self.initListenThread()
        printCmd('creating iw')
        self.size = float(size) if not size == 'default' else self.xdim
        
        self.iw = ImageWindow(image=self.projection, dim=self.size, label="Projection")
        self.iw.root.protocol("WM_DELETE_WINDOW", self.exitClient)
        self.iw.root.mainloop()
                

    def rotate(self, rot, tilt, psi):
        printCmd('image.projectVolumeDouble')
        self.fourierprojector.projectVolume(self.projection, rot, tilt, psi)
        printCmd('flipud')
        self.vol = flipud(getImageData(self.projection))
        printCmd('iw.updateData')
        self.iw.updateData(self.vol)
        printCmd('end rotate')
    
        
    def exit(self):
        XmippChimeraClient.exit(self)
        if not (self.iw is None):
            self.iw.root.destroy()
            
            
    def answer(self, msg):
        XmippChimeraClient.answer(self, msg)
        if msg == 'motion_stop':
            data = loads(self.client.recv())#wait for data
            printCmd('reading motion')
            self.motion = array(data)
            printCmd('getting euler angles')
            rot, tilt, psi = Euler_matrix2angles(self.motion)
            printCmd('calling rotate')  
            self.rotate(rot, tilt, psi)
            
            
    def exitClient(self):
        self.client.send('exit_client')
        self.exit()
        
        
     
            
def printCmd(cmd):
        pass
        #timeformat = "%S.%f" 
        #print datetime.now().strftime(timeformat) + ' %s'%cmd

