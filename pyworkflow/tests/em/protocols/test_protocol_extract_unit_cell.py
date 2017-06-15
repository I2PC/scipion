# ***************************************************************************
# * Authors:    Roberto Marabini (roberto@cnb.csic.es)
# *
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
# ***************************************************************************/
import os
import numpy as np
import math
import collections
from tempfile import mkstemp
from pyworkflow.utils import runJob
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.packages.xmipp3 import getEnviron
from pyworkflow.em.packages.xmipp3.constants import XMIPP_SYM_NAME
from pyworkflow.em.protocol import ProtImportVolumes
from pyworkflow.em.packages.xmipp3.protocol_extract_unit_cell import XmippProtExtractUnit
from pyworkflow.em.constants import SYM_I222r, SYM_I222, SCIPION_SYM_NAME, SYM_In25, SYM_In25r,\
     SYM_CYCLIC
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms_base import NMA_MASK_THRE
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms import XmippProtConvertToPseudoAtoms
from icosahedron import *
from pyworkflow.em.packages.ccp4.convert import Ccp4Header
from pyworkflow.em.data import Transform

OFFSET = 22.5

def generate_ico(sym, mode, f, ):
    icosahedron = Icosahedron(orientation=sym)
    center = []
    x = 0.
    y = 0.
    z = 0.
    if mode == 'xmipp':
        f.write("sph  + 1. %.3f %.3f %.3f .65\n" % (x, y, z))
    else:
        f.write('.sphere %.3f %.3f %.3f .65\n' % (x, y, z))
    # print 5fold points
    for vertice in icosahedron.get_vertices():
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (vertice[0], vertice[1], vertice[2]))
        else:
            f.write('.sphere %.3f %.3f %.3f .15\n' % (vertice[0], vertice[1], vertice[2]))

    # print 3fold points
    if mode == 'xmipp':
        pass
    else:
        f.write(""".comment 3fold
    .color yellow
    """)

    for _3fold in icosahedron.get_3foldAxis():
        x, y, z = _3fold
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .10\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .10\n' % (x, y, z))

    # print 2fold points
    if mode == 'xmipp':
        pass
    else:
        f.write(""".comment 2fold
    .color green
    """)
    for _2fold in icosahedron.get_2foldAxis():
        x, y, z = _2fold

        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .09\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .09\n' % (x, y, z))


def generate_cyclic(order, offset, mode, f, ):
    center = []
    z_value = [-.45, 0]
    if mode == 'xmipp':
        f.write("cyl  +  1. 0 0 +.45 1.0 1.0 .2 0 0 0\n")
        f.write("cyl  + -1. 0 0 +.45 0.8 0.8 .2 0 0 0\n")
        f.write("cyl  +  1. 0 0 +.45 0.5 0.5 .2 0 0 0\n")
        f.write("cyl  + -1. 0 0 +.45 0.5 0.5 .2 0 0 0\n")
    else:
        f.write('.cylinder 0 0 +.46 0 0 +.44 1 open\n')
        f.write('.cylinder 0 0 +.46 0 0 +.44 .5 open\n')
    for z in z_value:
        x = 0.
        y = 0.
        if mode == 'xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (x, y, z))
        else:
            f.write('.sphere %.3f %.3f %.3f .15\n' % (x, y, z))
        for point in range(order):
            x= math.cos(2*point*math.pi/order + offset)
            y= math.sin(2*point*math.pi/order + offset)
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .15\n" % (x,y,z))
            else:
                f.write('.sphere %.3f %.3f %.3f .15\n' % (x,y,z))
        for point in range(order):
            x= math.cos(2*point*math.pi/order + offset)/2
            y= math.sin(2*point*math.pi/order + offset)/2
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .10\n" % (x,y,z))
            else:
                f.write('.sphere %.3f %.3f %.3f .10\n' % (x,y,z))
        for point in range(order):
            x= math.cos(2*point*math.pi/order + offset)/4
            y= math.sin(2*point*math.pi/order + offset)/4
            if mode == 'xmipp':
                f.write("sph  + 1. %.3f %.3f %.3f .05\n" % (x,y,z))
            else:
                f.write('.sphere %.3f %.3f %.3f .05\n' % (x,y,z))



def generate(sym='I2n5', mode='xmipp',suffix="_i2", offset=0.):
    offset = math.radians(offset)
    symPreffix = sym[:1]
    symSuffix  = sym[1:]
    if mode=='xmipp':
        suffix += ".feat"
    else:
        suffix += ".bild"

    (fd, filename) = mkstemp(suffix=suffix)
    f = os.fdopen(fd, "w")


    if mode=='xmipp':
        f.write("""# Phantom description file, (generated with phantom help)
# General Volume Parameters:
#      Xdim      Ydim      Zdim   Background_Density Scale
       3 3 3 0 60
# Feature Parameters:
#Type  +/=  Density X_Center Y_Center Z_Center
""")
    else:
        f.write(""".scale 60
.comment system of coordinates
.color 1 0 0       
.arrow 0 0 0 2 0 0 0.025
.color 0 1 0
.arrow 0 0 0 0 2 0 0.025
.color 0 0 1
.arrow 0 0 0 0 0 2 0.025
.comment vertices
.color red
""")
    if symPreffix == 'I':
        generate_ico(symSuffix, mode, f)
    elif symPreffix == 'C':
        generate_cyclic(int(symSuffix), offset,  mode, f)

    f.close()
    return filename

class TestProtModelBuilding(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.mode = 'xmipp'
        #Cyclic
        cls.symOrder=8#phantom symmetry order for CN and similars

        #general
        cls.innerRadius = None
        cls.outerRadius= None
        setupTestProject(cls)
        cls.filename = {}
        cls.box = {}

    def test_extractunitCellcyclic(self):
        self.innerRadius = 0.
        self.outerRadius = 40.
        self.filename[SYM_CYCLIC] = generate(SCIPION_SYM_NAME[SYM_CYCLIC][:1]+str(self.symOrder) ,
                                             'xmipp', XMIPP_SYM_NAME[SYM_CYCLIC][:1]+str(self.symOrder))
        self.box[SYM_CYCLIC] = (76,49,81)
        self.extractunitCell(SYM_CYCLIC)
        self.filename[SYM_CYCLIC] = generate(SCIPION_SYM_NAME[SYM_CYCLIC][:1]+str(self.symOrder) ,
                                             'xmipp', XMIPP_SYM_NAME[SYM_CYCLIC][:1]+str(self.symOrder),OFFSET)
        self.box[SYM_CYCLIC] = (68,68,81)
        self.extractunitCell(SYM_CYCLIC, OFFSET)

    def test_extractunitCellIco(self):
        self.innerRadius = 37.
        self.outerRadius = 79.
        self.filename[SYM_I222] = generate(SCIPION_SYM_NAME[SYM_I222] , 'xmipp', XMIPP_SYM_NAME[SYM_I222])
        self.box[SYM_I222] = (93,91,53)
        self.filename[SYM_I222r] = generate(SCIPION_SYM_NAME[SYM_I222r] , 'xmipp', XMIPP_SYM_NAME[SYM_I222r])
        self.box[SYM_I222r] = (91,68,53)
        self.filename[SYM_In25] = generate(SCIPION_SYM_NAME[SYM_In25], 'xmipp', XMIPP_SYM_NAME[SYM_In25])
        self.box[SYM_In25] = (67,93,116)
        self.filename[SYM_In25r] = generate(SCIPION_SYM_NAME[SYM_In25r], 'xmipp', XMIPP_SYM_NAME[SYM_In25r])
        self.box[SYM_In25r] = (67,68,54)

        self.extractunitCell(SYM_I222)#no crowther 222
        self.extractunitCell(SYM_I222r)#crowther 222
        self.extractunitCell(SYM_In25)
        self.extractunitCell(SYM_In25r)

    def extractunitCell(self, sym, offset=0):

        """ extract unit cell from icosahedral pahntom
            using xmipp_i2 symmetry
        """
        # create phantom (3D map)
        _samplingRate = 1.34
        _, outputFile = mkstemp(suffix=".mrc")
        command = "xmipp_phantom_create "
        args    = " -i %s"% self.filename[sym]
        args += " -o %s"%outputFile
        runJob(None, command, args,env=getEnviron())
        ccp4header = Ccp4Header(outputFile, readHeader= True)
        ccp4header.setSampling(_samplingRate)
        x,y,z = ccp4header.getDims()
        t=Transform()
        t.setShifts(-x/2, -y/2, -z/2)
        ccp4header.setOffset(t)
        ccp4header.writeHeader()

        #import volume
        args = {'filesPath': outputFile,
                'filesPattern': '',
                'samplingRate': _samplingRate,
                'copyFiles': True,
                }
        prot = self.newProtocol(ProtImportVolumes, **args)
        prot.setObjLabel('import volume(%s)'% XMIPP_SYM_NAME[sym])
        self.launchProtocol(prot)

        # execute protocol extract unitCell
        args = {'inputVolumes': prot.outputVolume,
                'symmetryGroup': sym,
                'symmetryOrder': self.symOrder,
                'innerRadius': self.innerRadius,
                'outerRadius': self.outerRadius,
                'expandFactor': .2,
                'offset': offset
                }
        prot = self.newProtocol(XmippProtExtractUnit, **args)
        prot.setObjLabel('extract unit cell')
        self.launchProtocol(prot)

        #check results
        ih = ImageHandler()
        xdim, ydim, zdim, ndim = ih.getDimensions(prot.outputVolume.getFileName())
        self.assertEqual(xdim,self.box[sym][0])
        self.assertEqual(ydim,self.box[sym][1])
        self.assertEqual(zdim,self.box[sym][2])

        #create pdb file
        args = {'inputStructure':prot.outputVolume,
                'maskMode':NMA_MASK_THRE,
                'maskThreshold':0.5,
                'pseudoAtomRadius':1.5
                }
        prot = self.newProtocol(XmippProtConvertToPseudoAtoms, **args)
        prot.setObjLabel('get pdb')
        self.launchProtocol(prot)

        #check results
        filenamePdb = prot._getPath('pseudoatoms.pdb')
        self.assertTrue(os.path.isfile(filenamePdb))
        #delete temporary files
        os.remove(self.filename[sym])
        os.remove(outputFile)


