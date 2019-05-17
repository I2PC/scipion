# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from __future__ import print_function

import os
import ast
import shlex
import socket
import sys
from multiprocessing.connection import Client
from threading import Thread
from numpy import flipud

import pyworkflow as pw
import pyworkflow.utils as pwutils
from pyworkflow.viewer import CommandView, Viewer, DESKTOP_TKINTER
from pyworkflow.gui.matplotlib_image import ImageWindow
from pyworkflow.em.constants import (
    SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL, SYM_OCTAHEDRAL, SYM_I222,
    SYM_I222r, SYM_In25, SYM_In25r, SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r)
import pyworkflow.em.metadata as md
from pyworkflow.em.data import AtomStruct, PdbFile
from pyworkflow.em.convert import ImageHandler


import showj
import xmippLib


chimeraPdbTemplateFileName = "chimeraOut%04d.pdb"
chimeraMapTemplateFileName = "chimeraOut%04d.mrc"
chimeraScriptFileName = "chimeraScript.py"
sessionFile = "SESSION.py"

symMapperScipionchimera = {}
symMapperScipionchimera[SYM_CYCLIC] = "Cn"
symMapperScipionchimera[SYM_DIHEDRAL] = "Dn"
symMapperScipionchimera[SYM_TETRAHEDRAL] = "T"
symMapperScipionchimera[SYM_OCTAHEDRAL] = "O"
symMapperScipionchimera[SYM_I222] = "222"
symMapperScipionchimera[SYM_I222r] = "222r"
symMapperScipionchimera[SYM_In25] = "n25"
symMapperScipionchimera[SYM_In25r] = "n25r"
symMapperScipionchimera[SYM_I2n3] = "2n3"
symMapperScipionchimera[SYM_I2n3r] = "2n3r"
symMapperScipionchimera[SYM_I2n5] = "2n5"
symMapperScipionchimera[SYM_I2n5r] = "2n5r"


class Chimera:
    """ Helper class to execute chimera and handle its environment. """

    # Map symmetries from Scipion convention to Chimera convention
    _symmetryMap = {
        SYM_CYCLIC: 'Cn',
        SYM_DIHEDRAL: 'Dn',
        SYM_TETRAHEDRAL: 'T',
        SYM_OCTAHEDRAL: 'O',
        SYM_I222: '222',
        SYM_I222r: '222r',
        SYM_In25: 'n25',
        SYM_In25r: 'n25r',
        SYM_I2n3: '2n3',
        SYM_I2n3r: '2n3r',
        SYM_I2n5: '2n5',
        SYM_I2n5r: '2n5r'
    }

    @classmethod
    def getSymmetry(cls, scipionSym):
        """ Return the equivalent Chimera symmetry from Scipion one. """
        return cls._symmetryMap[scipionSym]

    @classmethod
    def getHome(cls):
        return os.environ.get('CHIMERA_HOME', None)

    @classmethod
    def getEnviron(cls):
        """ Return the proper environ to launch chimera.
        CHIMERA_HOME variable is read from the ~/.config/scipion.conf file.
        """
        environ = pwutils.Environ(os.environ)
        environ.set('PATH', os.path.join(cls.getHome(), 'bin'),
                    position=pwutils.Environ.BEGIN)

        if "REMOTE_MESA_LIB" in os.environ:
            environ.set('LD_LIBRARY_PATH', os.environ['REMOTE_MESA_LIB'],
                        position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def getProgram(cls, progName="chimera"):
        """ Return the program binary that will be used. """
        home = cls.getHome()
        if home is None:
            return None
        return os.path.join(home, 'bin', os.path.basename(progName))

    @classmethod
    def runProgram(cls, program=None, args=""):
        """ Internal shortcut function to launch chimera program. """
        prog = program or cls.getProgram()
        pwutils.runJob(None, prog, args, env=cls.getEnviron())

    @classmethod
    def createCoordinateAxisFile(cls, dim, bildFileName="/tmp/axis.bild",
                                 sampling=1, r1=0.1):
        """ Create a coordinate system, Along each dimension
        we place a small sphere in the negative axis. In this way
        chimera shows the system of coordinates origin in the
        window center"""

        ff = open(bildFileName, "w")
        arrowDict = {}
        arrowDict["x"] = arrowDict["y"] = arrowDict["z"] = \
            sampling * dim * 3. / 4.
        arrowDict["r1"] = r1 * dim / 50.
        arrowDict["r2"] = 4 * r1
        arrowDict["rho"] = 0.75  # axis thickness

        ff.write(".color 1 0 0\n"
                 ".arrow 0 0 0 %(x)f 0 0 %(r1)f %(r2)f %(rho)f\n"
                 ".color 0 0 0\n.sphere -%(x)f 0 0 0.00001\n"
                 ".color 1 1 0\n"
                 ".arrow 0 0 0 0 %(y)f 0 %(r1)f %(r2)f %(rho)f\n"
                 ".color 0 0 0\n.sphere 0 -%(y)f 0 0.00001\n"
                 ".color 0 0 1\n"
                 ".arrow 0 0 0 0 0 %(z)f %(r1)f %(r2)f %(rho)f\n"
                 ".color 0 0 0\n.sphere 0 0 -%(z)f 0.00001\n" %
                 arrowDict)
        ff.close()


def printCmd(cmd):
    """ Debug funtion. """
    pass
    # print cmd



class ChimeraClient:

    def openVolumeOnServer(self, volume, sendEnd=True):
        self.send('open_volume', volume)
        if self.voxelSize is not None:
            self.send('voxel_size', self.voxelSize)
        if sendEnd:
            self.client.send('end')

    def __init__(self, volfile, sendEnd=True, **kwargs):
        if volfile.endswith('.mrc'):
            volfile += ':mrc'

        self.kwargs = kwargs
        if volfile is None:
            raise ValueError(volfile)
        if '@' in volfile:
            [index, file] = volfile.split('@')
        else:
            file = volfile
        if ':' in file:
            file = file[0: file.rfind(':')]
        if not os.path.exists(file):
            raise Exception("File %s does not exists" % file)

        self.volfile = volfile
        self.voxelSize = self.kwargs.get('voxelSize', None)
        # ChimeraServer is the basic server. There are other
        # than inherit from it.
        serverName = self.kwargs.get('ChimeraServer', 'ChimeraServer')
        self.address = ''
        self.port = pwutils.getFreePort()

        serverfile = pw.join('em', 'viewers', 'chimera_server.py')
        command = CommandView("chimera --script '%s %s %s' &" %
                              (serverfile, self.port, serverName),
                              env=Chimera.getEnviron(),).show()
        self.authkey = 'test'
        self.client = Client((self.address, self.port), authkey=self.authkey)
        self.initVolumeData()
        # self.openVolumeOnServer(self.vol,sendEnd)
        self.openVolumeOnServer(self.vol)
        self.initListenThread()

    def send(self, cmd, data):
        self.client.send(cmd)
        self.client.send(data)

    def initListenThread(self):
            self.listen_thread = Thread(name="ChimeraCli.listenTh",
                                        target=self.listen)
            # self.listen_thread.daemon = True
            self.listen_thread.start()

    def listen(self):

        self.listen = True
        try:
            while self.listen:
                msg = self.client.recv()
                self.answer(msg)
        except EOFError:
            print ('Lost connection to server')
        finally:
            self.exit()

    def exit(self):
        self.client.close()

    def initVolumeData(self):
        self.image = xmippLib.Image(self.volfile)
        self.image.convert2DataType(md.DT_DOUBLE)
        self.xdim, self.ydim, self.zdim, self.n = self.image.getDimensions()
        self.vol = self.image.getData()

    def answer(self, msg):
        if msg == 'exit_server':
            self.listen = False


class ChimeraAngDistClient(ChimeraClient):

    def __init__(self, volfile, **kwargs):
        self.angularDistFile = kwargs.get('angularDistFile', None)

        if self.angularDistFile:
            fileName = self.angularDistFile
            if '@' in self.angularDistFile:
                fileName = self.angularDistFile.split("@")[1]
            if not (os.path.exists(fileName)):  # check blockname:
                raise Exception("Path %s does not exists" %
                                self.angularDistFile)
        self.spheresColor = kwargs.get('spheresColor', 'red')
        spheresDistance = kwargs.get('spheresDistance', None)
        spheresMaxRadius = kwargs.get('spheresMaxRadius', None)
        ChimeraClient.__init__(self, volfile, sendEnd=False, **kwargs)
        self.spheresDistance = float(spheresDistance) \
            if spheresDistance else 0.75 * max(self.xdim, self.ydim, self.zdim)
        self.spheresMaxRadius = float(spheresMaxRadius) \
            if spheresMaxRadius else 0.02 * self.spheresDistance
        self.loadAngularDist(True)

    def loadAngularDist(self,  sendEnd=True):
        if self.angularDistFile:
            self.readAngularDistFile()
            self.send('command_list', self.angulardist)
        if sendEnd:
            self.client.send('end')

    def readAngularDistFile(self):
        angleRotLabel = md.MDL_ANGLE_ROT
        angleTiltLabel = md.MDL_ANGLE_TILT
        anglePsiLabel = md.MDL_ANGLE_PSI
        mdAngDist = md.MetaData(self.angularDistFile)
        if not mdAngDist.containsLabel(md.MDL_ANGLE_PSI):
            anglePsiLabel = None
            if mdAngDist.containsLabel(md.RLN_ORIENT_PSI):
                angleRotLabel = md.RLN_ORIENT_ROT
                angleTiltLabel = md.RLN_ORIENT_TILT
                anglePsiLabel = md.RLN_ORIENT_PSI

        if not mdAngDist.containsLabel(md.MDL_WEIGHT):
            mdAngDist.fillConstant(md.MDL_WEIGHT, 1.)

        maxweight = mdAngDist.aggregateSingle(md.AGGR_MAX, md.MDL_WEIGHT)
        minweight = mdAngDist.aggregateSingle(md.AGGR_MIN, md.MDL_WEIGHT)
        interval = maxweight - minweight

        self.angulardist = []
        x2 = self.xdim/2
        y2 = self.ydim/2
        z2 = self.zdim/2
        # cofr does not seem to work!
        # self.angulardist.append('cofr %d,%d,%d'%(x2,y2,z2))
        for id in mdAngDist:
            rot = mdAngDist.getValue(angleRotLabel, id)
            tilt = mdAngDist.getValue(angleTiltLabel, id)
            psi = mdAngDist.getValue(anglePsiLabel, id) if anglePsiLabel else 0
            weight = mdAngDist.getValue(md.MDL_WEIGHT, id)
            # Avoid zero division
            weight = 0 if interval == 0 else (weight - minweight)/interval
            weight = weight + 0.5  # add 0.5 to avoid cero weight
            x, y, z = xmippLib.Euler_direction(rot, tilt, psi)
            radius = weight * self.spheresMaxRadius

            x = x * self.spheresDistance + x2
            y = y * self.spheresDistance + y2
            z = z * self.spheresDistance + z2
            command = 'shape sphere radius %s center %s,%s,%s color %s ' % \
                      (radius, x, y, z, self.spheresColor)
            self.angulardist.append(command)
            # printCmd(command)


class ChimeraVirusClient(ChimeraClient):

    def __init__(self, volfile, **kwargs):
        self.h = kwargs.get('h', 5)
        if self.h is None:
                self.h = 5
        self.k = kwargs.get('k', 0)
        if self.k is None:
                self.k = 0
        self.sym = kwargs.get('sym', 'i222r')
        if self.sym is None:
                self.sym = 'i222r'
        self.radius = kwargs.get('radius', 100.)
        if self.radius is None:
                self.radius = 100.
        self.spheRadius = kwargs.get('spheRadius', 1.5)
        if self.spheRadius is None:
                self.spheRadius = 1.5
        self.color = kwargs.get('color', 'red')
        if self.color is None:
                self.color = 'red'
        self.linewidth = kwargs.get('linewidth', 1)
        if self.linewidth is None:
                self.linewidth = 1
        self.sphere = kwargs.get('sphere', 0)
        if self.sphere is None:
                self.sphere = 0
        self.shellRadius = kwargs.get('shellRadius', self.spheRadius)
        if self.shellRadius is None:
                self.shellRadius = self.spheRadius
        kwargs['ChimeraServer'] = 'ChimeraVirusServer'
        ChimeraClient.__init__(self, volfile, **kwargs)

    def hkcageCommand(self):
        command = 'hkcage '
        command += '%d %d ' % (self.h, self.k)
        command += 'radius  %f ' % self.radius
        command += 'orientation %s ' % self.sym
        command += 'replace true '
        command += 'color %s ' % self.color
        command += 'linewidth %d ' % self.linewidth
        command += 'sphere %f ' % self.sphere
        return command

    def openVolumeOnServer(self, volume):
        ChimeraClient.openVolumeOnServer(self, volume, sendEnd=False)
        commandList = [self.hkcageCommand()]
        self.send('hk_icosahedron_lattice', (self.h,
                                             self.k,
                                             self.radius,
                                             self.shellRadius,
                                             self.spheRadius,
                                             self.sym,
                                             self.sphere,
                                             self.color))
        # get here the list of vertexes, info will be pass by id later
        self.send('command_list', commandList)
        # endhere
        self.client.send('end')
        return

        # get va with sphere centers

        # va coordinates of  vertex of the triangles inside de canonical
        # triangle
        msg1 = self.client.recv()
        self.va = self.client.recv()
        # self.listToBild(self.va,1.6,msg1+'.bild')
        self.client.send('end')

    # this is an auxiliary file that should be deleted
    def listToBild(self, points, radius, file, color='0 0 1'):
        f = open(file, 'w')
        for point in points:
            print ("\n.color", color, "\n.sphere", point[0]+128.,
                   point[1]+128., point[2]+128., radius, file=f)

    def answer(self, msg):
        ChimeraClient.answer(self, msg)
        if msg == 'motion_stop':
            data = self.client.recv()  # wait for data
            printCmd('reading motion')
            self.motion = data
            printCmd('getting euler angles')
            rot, tilt, psi = xmippLib.Euler_matrix2angles(self.motion)
            printCmd('calling rotate')
            self.rotate(rot, tilt, psi)
        elif msg == 'id':
            id = self.client.recv()
            self.listToBild([self.va[id-1]], 2.5, 'id.bild')
        else:
            pass

    def listenShowJ(self):
        """ This function with a very confusing name send messages to
        the chimera server"""

        while True:
            try:
                (clientsocket, address) = self.serversocket.accept()
                # Should be a single message, so no loop
                msg = clientsocket.recv(1024)
                tokens = shlex.split(msg)
                cmd = tokens[0]
                if cmd == 'rotate':
                    rot = float(tokens[1])
                    tilt = float(tokens[2])
                    psi = float(tokens[3])

                    matrix = xmippLib.Euler_angles2matrix(rot, tilt, psi)
                elif cmd == 'rotate_matrix':
                    matrixString = tokens[1]
                    matrix = ast.literal_eval(matrixString)
                    matrix = [matrix[0][:3], matrix[1][:3], matrix[2][:3]]
                self.client.send('rotate')
                self.client.send(matrix)
                clientsocket.close()

            except EOFError:
                print ('Lost connection to client')


class ChimeraProjectionClient(ChimeraAngDistClient):

    def __init__(self, volfile, **kwargs):
        print("on chimera projection client")
        ChimeraAngDistClient.__init__(self, volfile, **kwargs)
        self.projection = xmippLib.Image()
        self.projection.setDataType(md.DT_DOUBLE)
        # 0.5 ->  Niquiest frequency
        # 2 -> bspline interpolation
        size = self.kwargs.get('size', None)
        defaultSize = self.xdim if self.xdim > 128 else 128
        self.size = size if size else defaultSize
        paddingFactor = self.kwargs.get('paddingFactor', 1)
        maxFreq = self.kwargs.get('maxFreq', 0.5)
        splineDegree = self.kwargs.get('splineDegree', 3)
        self.fourierprojector = xmippLib.FourierProjector(
            self.image, paddingFactor, maxFreq, splineDegree)
        self.fourierprojector.projectVolume(self.projection, 0, 0, 0)
        self.showjPort = self.kwargs.get('showjPort', None)
        self.iw = ImageWindow(filename=os.path.basename(volfile),
                              image=self.projection,
                              dim=self.size, label="Projection")
        self.iw.updateData(flipud(self.projection.getData()))
        if self.showjPort:
            self.showjThread = Thread(name="ChimeraProjClient",
                                      target=self.listenShowJ)
            self.showjThread.daemon = True
            self.showjThread.start()
        self.iw.root.protocol("WM_DELETE_WINDOW", self.exitClient)
        self.iw.show()

    def rotate(self, rot, tilt, psi):
        self.fourierprojector.projectVolume(self.projection, rot, tilt, psi)
        self.projectionData = flipud(self.projection.getData())

        if hasattr(self, 'iw'):
            # sometimes is not created and rotate is called
            self.iw.updateData(self.projectionData)

    def exit(self):
        ChimeraClient.exit(self)
        if hasattr(self, "iw"):
            self.iw.root.destroy()

    def answer(self, msg):
        ChimeraClient.answer(self, msg)
        if msg == 'motion_stop':
            data = self.client.recv()  # wait for data
            printCmd('reading motion')
            self.motion = data
            printCmd('getting euler angles')
            rot, tilt, psi = xmippLib.Euler_matrix2angles(self.motion)
            printCmd('calling rotate')
            self.rotate(rot, tilt, psi)

    def exitClient(self):  # close window before volume loaded
        if not self.listen:
            sys.exit(0)

    def initListenThread(self):
        self.listen_thread = Thread(name="ChimeraProjectionCli.initListen",
                                    target=self.listen)
        self.listen_thread.daemon = True
        self.listen_thread.start()

    def listenShowJ(self):
        self.serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.serversocket.bind(('', self.showjPort))
        # become a server socket
        self.serversocket.listen(1)

        while True:
            try:
                (clientsocket, address) = self.serversocket.accept()
                # should be a single message, so no loop
                msg = clientsocket.recv(1024)
                tokens = shlex.split(msg)
                cmd = tokens[0]
                if cmd == 'rotate':
                    rot = float(tokens[1])
                    tilt = float(tokens[2])
                    psi = float(tokens[3])

                    matrix = xmippLib.Euler_angles2matrix(rot, tilt, psi)
                elif cmd == 'rotate_matrix':
                    matrixString = tokens[1]
                    matrix = ast.literal_eval(matrixString)
                    matrix = [matrix[0][:3], matrix[1][:3], matrix[2][:3]]
                self.client.send('rotate')
                self.client.send(matrix)
                clientsocket.close()

            except EOFError:
                print ('Lost connection to client')


class ChimeraView(CommandView):
    """ View for calling an external command. """
    def __init__(self, inputFile, **kwargs):
        CommandView.__init__(self, 'chimera "%s" &' % inputFile,
                             env=Chimera.getEnviron(), **kwargs)


class ChimeraClientView(CommandView):
    """ View for calling an external command. """
    def __init__(self, inputFile, **kwargs):
        self._inputFile = inputFile
        self._kwargs = kwargs

    def show(self):
        if self._kwargs.get('showProjection', False):
            ChimeraProjectionClient(self._inputFile, **self._kwargs)
        else:
            ChimeraAngDistClient(self._inputFile, **self._kwargs)


class ChimeraDataView(ChimeraClientView):

    def __init__(self, dataview, vol, viewParams={}, **kwargs):
        self.dataview = dataview
        self.showjPort = pwutils.getFreePort()
        self.dataview._viewParams[showj.CHIMERA_PORT] = self.showjPort
        self.dataview._viewParams[showj.MODE] = showj.MODE_MD
        self.dataview._viewParams[showj.INVERTY] = ''
        ChimeraClientView.__init__(self, vol.getFileName(),
                                   showProjection=True,
                                   showjPort=self.showjPort,
                                   voxelSize=vol.getSamplingRate())

    def show(self):
        self.dataview.show()
        ChimeraClientView.show(self)


class ChimeraViewer(Viewer):
    """ Wrapper to visualize PDB object with Chimera. """
    _environments = [DESKTOP_TKINTER]
    _targets = [AtomStruct, PdbFile]

    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    def visualize(self, obj, **kwargs):
        cls = type(obj)
        if issubclass(cls, AtomStruct):
            # if attribute _chimeraScript exists then protocol
            # has create a script file USE IT
            if hasattr(obj, '_chimeraScript'):
                fn = obj._chimeraScript.get()
                ChimeraView(fn).show()
                return
            # if not create a script file with: coordinates axis, PDB and
            # volume (if available)
            else:
                fn = obj.getFileName()
                # check if tmp dir exists, if not use /tmp
                # tmp does not exists if you try to visualize something  (eye)
                # before irunning the protocol
                tmpPath=self.protocol._getTmpPath()
                if not os.path.exists(tmpPath):
                    tmpPath = "/tmp"
                fnCmd = os.path.join(tmpPath, "chimera.cmd")
                f = open(fnCmd, 'w')
                f.write("cofr 0,0,0\n")  # set center of coordinates
                if obj.hasVolume():
                    volID = 0
                    volumeObject = obj.getVolume()
                    dim = volumeObject.getDim()[0]
                    sampling = volumeObject.getSamplingRate()
                    f.write("open %s\n" % os.path.abspath(
                        ImageHandler.removeFileType(volumeObject.getFileName())))
                    f.write("volume #%d style surface voxelSize %f\n"
                            % (volID, sampling))
                    x, y, z = volumeObject.getShiftsFromOrigin()
                    f.write("volume #%d origin %0.2f,%0.2f,%0.2f\n"
                            % (volID, x, y, z))
                else:
                    dim = 150  # eventually we will create a PDB library that
                               # computes PDB dim
                    sampling = 1.
                # Construct the coordinate file
                bildFileName = os.path.abspath(
                    os.path.join(tmpPath, "axis.bild"))
                Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=bildFileName,
                                         sampling=sampling)
                f.write("open %s\n" % bildFileName)
                f.write("open %s\n" % os.path.abspath(fn))
                f.close()
                ChimeraView(fnCmd).show()
            # FIXME: there is an asymmetry between ProtocolViewer and Viewer
            # for the first, the visualize method return a list of View's
            # (that are shown)
            # for the second, the visualize directly shows the objects.
            # the first approach is better
        else:
            raise Exception('ChimeraViewer.visualize: '
                            'can not visualize class: %s' % obj.getClassName())
