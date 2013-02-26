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

import xmipp
from multiprocessing.connection import Client
from protlib_xmipp import getImageData
from protlib_gui_figure import ImageWindow
from pickle import dumps, loads
from threading import Thread
from os import system
from numpy import array, ndarray, flipud

class ChimeraClient:
    
    def __init__(self, projexplorer):
        self.projexplorer = projexplorer
        #print 'volfile: ' + self.volfile
        self.address = ''
        self.port = 6000
        self.authkey = 'test'
        self.client = Client((self.address, self.port), authkey=self.authkey)
        
        self.listen_thread = Thread(target=self.listen)
        self.listen_thread.daemon = True
        self.listen_thread.start()
        

        
    def openVolume(self, Z):
        self.client.send(Z)
        
         
    def listen(self):
        
        try:
            while True:
                #print 'on client loop'
                msg = self.client.recv()
                cmd = loads(msg) 
                if cmd == 'exit':
                    print 'exit msg received'
                    break
                else:
                    self.motion = array(cmd)
                    rot1, tilt1, psi1 = xmipp.Euler_matrix2angles(self.motion)  
                    self.projexplorer.rotate(rot1, tilt1, psi1)
        except EOFError:
            print 'Lost connection to server'
        finally:
            #print 'closing client'
            self.client.close()#close connection
            self.projexplorer.close()

    



class XmippProjectionExplorer:
    
    def __init__(self, volfile):
        self.volfile = volfile
        self.initProjection()
        self.client = ChimeraClient(self)
        self.client.openVolume(self.getVolumeData())
        self.iw = ImageWindow(image=self.projection, label="Projection")
        self.iw.root.mainloop()
     
        
    def initProjection(self):
        self.image = xmipp.Image(self.volfile)
        self.image.convert2DataType(xmipp.DT_DOUBLE)
        self.projection = self.image.projectVolumeDouble(0, 0, 0)
        
    def rotate(self, rot, tilt, psi):
        self.projection = self.image.projectVolumeDouble(rot, tilt, psi)
        Z = flipud(getImageData(self.projection))
        self.iw.updateData(Z)
        
    def getVolumeData(self):
        xdim, ydim, zdim, n = self.image.getDimensions()
        Z = getImageData(self.image)
        return Z

    def close(self):
        self.iw.root.destroy()
