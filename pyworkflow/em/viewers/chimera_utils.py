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
"""
This module contains converter functions that will serve to:
1. define ccp4 environ
TODO:
2. Read/Write CCP4 specific files
"""

import os
import pyworkflow.utils as pwutils
from pyworkflow.em.constants import SYM_CYCLIC, SYM_DIHEDRAL, \
    SYM_TETRAHEDRAL, SYM_OCTAHEDRAL, SYM_I222, SYM_I222r, SYM_In25, SYM_In25r
from pyworkflow.utils import  Environ

#
chimeraPdbTemplateFileName = "scipionOut%04d.pdb"
chimeraMapTemplateFileName = "scipionOut%04d.mrc"
chimeraScriptFileName = "chimeraScript.py"

symMapperScipionchimera = {}
symMapperScipionchimera[SYM_CYCLIC] = "Cn"
symMapperScipionchimera[SYM_DIHEDRAL] = "Dn"
symMapperScipionchimera[SYM_TETRAHEDRAL] = "T"
symMapperScipionchimera[SYM_OCTAHEDRAL] = "O"
symMapperScipionchimera[SYM_I222] = "222"
symMapperScipionchimera[SYM_I222r] = "222r"
symMapperScipionchimera[SYM_In25] = "n25"
symMapperScipionchimera[SYM_In25r] = "n25r"


def getEnviron():
    return getChimeraEnviron()


def runChimeraProgram(program, args=""):
    """ Internal shortcut function to launch chimera program. """
    env = getEnviron()
    pwutils.runJob(None, program, args, env=env)


def getProgram(progName="chimera"):
    """ Return the program binary that will be used. """
    if 'CHIMERA_HOME' not in os.environ:
        return None
    return os.path.join(os.environ['CHIMERA_HOME'], 'bin',
                        os.path.basename(progName))


def createCoordinateAxisFile(dim, bildFileName="/tmp/axis.bild",
                             sampling=1, r1=0.1):
    ff = open(bildFileName, "w")
    arrowDict = {}
    arrowDict["x"] = arrowDict["y"] = arrowDict["z"] = \
        sampling * dim * 3. / 4.
    arrowDict["r1"] = r1 * dim / 50.
    arrowDict["r2"] = 4 * r1
    arrowDict["rho"] = 0.75  # axis thickness

    ff.write(".color 1 0 0\n"
             ".arrow 0 0 0 %(x)d 0 0 %(r1)f %(r2)f %(rho)f\n"
             ".color 1 1 0\n"
             ".arrow 0 0 0 0 %(y)d 0 %(r1)f %(r2)f %(rho)f\n"
             ".color 0 0 1\n"
             ".arrow 0 0 0 0 0 %(z)d %(r1)f %(r2)f %(rho)f\n" %
             arrowDict)
    ff.close()

def adaptOriginFromCCP4ToChimera(origin):
    return tuple(-1.0*x for x in origin)


def getChimeraEnviron():
    """ Return the proper environ to launch chimera.
    CHIMERA_HOME variable is read from the ~/.config/scipion.conf file.
    """
    environ = Environ(os.environ)
    environ.set('PATH', os.path.join(os.environ['CHIMERA_HOME'], 'bin'),
                position=Environ.BEGIN)

    if "REMOTE_MESA_LIB" in os.environ:
        environ.set('LD_LIBRARY_PATH', os.environ['REMOTE_MESA_LIB'],
                    position=Environ.BEGIN)
    return environ


# define scipion_write command
chimeraScriptHeader = '''
import os
def beep(time):
   """I simply do not know how to create a portable beep sound.
      This system call seems to work pretty well if you have sox
      installed"""
   try:
      command = "play --no-show-progress -n synth %f sin 880"%time
      os.system(command)
   except:
      pass

"""load the script in chimera for example chimera --script ScipioChimeraExt/
ChimeraExtension.py
from chimera command line type one of the following three options
scipionwrite
scipionwrite model #n
scipionwrite model #n refmodel #p saverefmodel 0/1
"""
def cmd_scipionWrite(scipionWrite,args):
  from Midas.midas_text import doExtensionFunc
'''

chimeraScriptMain = '''
  def scipionWrite(model="#%(pdbID)d",refmodel="#%(_3DmapId)d",
     saverefmodel=0):
     #get model (pdb) id
     modelId = int(model[1:])# model to write
     refModelId = int(refmodel[1:])# coordenate system refers to this model

     # get actual models
     model    = chimera.openModels.list()[modelId]
     refModel = chimera.openModels.list()[refModelId]

     # Save the PDB relative to the volume coordinate system
     # TODO: check if this Will work if the reference is a PDB?
     from Midas import write
     write(model, refModel, "%(pdbFileTemplate)s")
     # alternative way to save  the pdb file using a command
     #run('write relative #1 #0 pdb_path')

     # Save the map if sampling rate has been changed
     if saverefmodel:
         from VolumeData import save_grid_data
         save_grid_data(refModel.data, "%(chimeraMapTemplateFileName)s")
     beep(0.1)

  doExtensionFunc(scipionWrite,args)

from Midas.midas_text import addCommand
addCommand('scipionwrite', cmd_scipionWrite, help="http://scipion.cnb.csic.es")
'''


def createScriptFile(pdbID, _3DmapId,
                     pdbFileTemplate, mapFileTemplate, f):
    f.write(chimeraScriptHeader)
    d = {}
    d['pdbID'] = pdbID
    d['_3DmapId'] = _3DmapId
    d['pdbFileTemplate'] = pdbFileTemplate % 1
    d['chimeraMapTemplateFileName'] = mapFileTemplate % 1

    f.write(chimeraScriptMain % d)
