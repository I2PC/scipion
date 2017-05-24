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
from tempfile import mkstemp
from pyworkflow.utils import runJob
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.packages.xmipp3 import getEnviron
from pyworkflow.em.protocol import ProtImportVolumes
from pyworkflow.em.packages.xmipp3.protocol_extract_unit_cell import XmippProtExtractUnit
from pyworkflow.em.constants import SYM_I222r
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms_base import NMA_MASK_THRE
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms import XmippProtConvertToPseudoAtoms

def generate(fiveFold, mode='xmipp',suffix=".feat"):

    (fd, filename) = mkstemp(suffix=suffix)
    f = os.fdopen(fd, "w")
    #f.write(geometricPhantom)
    #f.close()

    if mode=='xmipp':
        f.write("""# Phantom description file, (generated with phantom help)
# General Volume Parameters:
#      Xdim      Ydim      Zdim   Background_Density Scale
       5 5 5 0 60
# Feature Parameters:
#Type  +/=  Density X_Center Y_Center Z_Center
""")
    else:
        f.write(""".scale 60
.comment five fold
.color red
""")
    counter=0
    for vertex in fiveFold:
        if mode=='xmipp':
            f.write("sph  + 1. %.3f %.3f %.3f .15\n"%(vertex[0], vertex[1], vertex[2]) )
        else:
            f.write(".cmov   %.3f %.3f %.3f\n%d\n"% (vertex[0], vertex[1], vertex[2], counter))
            f.write('.sphere %.3f %.3f %.3f .15\n'% (vertex[0], vertex[1], vertex[2]))
        counter += 1

    ### vertex that define a triangle
    threeFold  = []
    #zero
    threeFold.append(np.array([0,2,5]))
    threeFold.append(np.array([0,4,2]))
    threeFold.append(np.array([0,5,10]))
    threeFold.append(np.array([0,8,4]))
    threeFold.append(np.array([0,10,8]))
    #one
    threeFold.append(np.array([1,3,6]))
    threeFold.append(np.array([1,6,8]))
    threeFold.append(np.array([1,7,3]))
    threeFold.append(np.array([1,8,10]))
    threeFold.append(np.array([1,10,7]))
    #two
    threeFold.append(np.array([2,4,9]))
    threeFold.append(np.array([2,9,11]))
    threeFold.append(np.array([2,11,5]))
    #three
    threeFold.append(np.array([3,6,9]))
    threeFold.append(np.array([3,9,11]))
    threeFold.append(np.array([3,11,7]))
    #four
    threeFold.append(np.array([4,6,9]))
    threeFold.append(np.array([4,8,6]))
    #five
    threeFold.append(np.array([5,7,11]))
    threeFold.append(np.array([5,10,7]))

    if mode == "xmipp":
        pass
    else:
        f.write(""".comment three fold
.color yellow
""")
    for cont in threeFold:
        temp = (fiveFold[cont[0]]+fiveFold[cont[1]]+fiveFold[cont[2]])/3.
        if mode == "xmipp":
            f.write("sph  + 1. %.3f %.3f %.3f .15\n"%(temp[0], temp[1], temp[2]))
        else:
            f.write(".sphere   %.3f %.3f %.3f .15\n"%(temp[0], temp[1], temp[2]))
    ###
    contFiveFold = []
    for fold in threeFold:
        contFiveFold.append(np.array([fold[0],fold[1]]))
        contFiveFold.append(np.array([fold[1],fold[2]]))
        contFiveFold.append(np.array([fold[2],fold[0]]))
    if mode == "xmipp":
        pass
    else:
        f.write(""".comment two fold
.color green
""")
    for cont in contFiveFold:
        temp = (fiveFold[cont[0]]+fiveFold[cont[1]])/2.
        if mode == "xmipp":
            f.write("sph  + 1. %.3f %.3f %.3f .15\n"%(temp[0], temp[1], temp[2]))
        else:
            f.write(".sphere %.3f %.3f %.3f .15\n"%(temp[0], temp[1], temp[2]))


    if mode == "bild":
        f.write(""".comment edges
.color blue
""")
        for cont in contFiveFold:
            temp = (fiveFold[cont[0]]+fiveFold[cont[1]])/2.
            f.write(".cylinder %.3f %.3f %.3f %.3f %.3f %.3f .03\n"%(fiveFold[cont[0]][0],
         fiveFold[cont[0]][1],
         fiveFold[cont[0]][2],
         fiveFold[cont[1]][0],
         fiveFold[cont[1]][1],
         fiveFold[cont[1]][2]))

    f.close()
    return filename

class TestProtModelBuilding(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        #i2
        fiveFold = []
        fiveFold.append(np.array([+1.000,    0,    1.618]))
        fiveFold.append(np.array([+1.000,    0,   -1.618]))
        fiveFold.append(np.array([-1.000,    0,    1.618]))
        fiveFold.append(np.array([-1.000,    0,   -1.618]))
        #
        fiveFold.append(np.array([ 0.,      1.618,  +1.000]))
        fiveFold.append(np.array([ 0.,      -1.618, +1.000]))
        fiveFold.append(np.array([ 0.,      1.618,  -1.000]))
        fiveFold.append(np.array([ 0.,      -1.618, -1.000]))
        #
        fiveFold.append(np.array([ 1.618,  +1.000, 0.]))
        fiveFold.append(np.array([-1.618,  +1.000, 0.]))
        fiveFold.append(np.array([ 1.618,  -1.000, 0.]))
        fiveFold.append(np.array([-1.618,  -1.000, 0.]))

        cls.filename = generate(fiveFold, 'xmipp',"_i2.feat")

    def test_extractunitCell(self):
        """ extract unit cell from icosahedral pahntom
            using xmipp_i2 symmetry
        """
        # create phantom (3D map)
        _, outputFile = mkstemp(suffix=".vol")
        command = "xmipp_phantom_create "
        args    = " -i %s"%self.filename
        args += " -o %s"%outputFile
        runJob(None, command, args,env=getEnviron())

        #import volume
        args = {'filesPath': outputFile,
                'filesPattern': '',
                'samplingRate': 1.34,
                'copyFiles': True,
                }
        prot = self.newProtocol(ProtImportVolumes, **args)
        prot.setObjLabel('import volume')
        self.launchProtocol(prot)

        # execute protocol extract unitCell
        args = {'inputVolumes': prot.outputVolume,
                'symmetryGroup': SYM_I222r,
                'innerRadius': 80,
                'outerRadius': 124,
                'expandFactor': .2,
                'offset': 0.
                }
        prot = self.newProtocol(XmippProtExtractUnit, **args)
        prot.setObjLabel('extract unit cell')
        self.launchProtocol(prot)

        #check results
        ih = ImageHandler()
        xdim, ydim, zdim, ndim = ih.getDimensions(prot.outputVolume.getFileName())
        self.assertEqual(xdim,111)
        self.assertEqual(ydim,107)
        self.assertEqual(zdim,69)

        #create pdb file
        args = {'inputStructure':prot.outputVolume,
                'maskMode':NMA_MASK_THRE,
                'maskThreshold':1.,
                'pseudoAtomRadius':1.5
                }
        prot = self.newProtocol(XmippProtConvertToPseudoAtoms, **args)
        prot.setObjLabel('get pdb')
        self.launchProtocol(prot)
        #check results
        filenamePdb = prot._getPath('pseudoatoms.pdb')
        self.assertTrue(os.path.isfile(filenamePdb))
        #delete temporary files
        os.remove(self.filename)
        os.remove(outputFile)


