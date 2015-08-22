# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin  (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
import sys
import json

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
import pyworkflow.em as em



def getPyTomPaths():
    """ Return the list of paths that need to be included for use PyTom. """
    PYTOM_HOME = os.environ['PYTOM_HOME']
    suffixes = ['pytomc', 'pytomc/swigModules', 
                'pytomc/sh_alignment/frm/swig',                 
                'pytomc/libs/libtomc/libs',
                'external/lib', 'external/lib/python/site-packages']
    suffixPaths = [os.path.join(PYTOM_HOME, suffix) for suffix in suffixes]
    
    return [os.path.dirname(PYTOM_HOME)] + suffixPaths


def addPyTomPaths():
    """ Add all path to current sys.path """
    for p in getPyTomPaths():
        sys.path.append(p)
    
    
def getEnviron():
    """ Create the needed environment for Pytom programs. """
    environ = pwutils.Environ(os.environ)
    if 'PYTOM_HOME' in os.environ:
        PYTOM_HOME = os.environ['PYTOM_HOME']
    else:
        print ('Warning: Cannot find environ variable PYTOM_HOME. Using '
               '$SCIPION_HOME/software/em/pytom instead.')
        PYTOM_HOME = '%s/software/em/pytom' % os.environ['SCIPION_HOME']
    environ.update({
#                    'LD_LIBRARY_PATH': os.pathsep.join([
#                                       os.path.join(PYTOM_HOME, 'pytomc/libs/libtomc/libs'),
#                                       os.path.join(PYTOM_HOME, 'pytomc/sh_alignment/SpharmonicKit27')]),
                    'PATH': os.path.join(PYTOM_HOME, 'bin'),
                    'PYTHONPATH': os.pathsep.join(getPyTomPaths())
                    }, 
                   position=pwutils.Environ.BEGIN)
    return environ


def writeSetOfVolumes(volSet, volXml, volDir):
    """ Convert a SetOfVolumes to a xml file used by PyTom.
    The volumes will be converted to .mrc format if not are '.em' or '.mrc' 
    Params:
        volSet: input SetOfVolumes.
        volXml: filename where to write the xml file.
        volDir: where to create links or copies (converted to mrc).
    """
    # Add to the path the root to pytom
    backupPath = list(sys.path)
    addPyTomPaths()
    
    from pytom.basic.structures import Particle, ParticleList, Wedge, SingleTiltWedge
    from pytom.score.score import Score, PeakPrior, xcfScore
    from pytom.frm.FRMAlignment import FRMScore
    
    w = SingleTiltWedge()
    #s = PeakPrior()
    pl = ParticleList()
    ih = em.convert.ImageHandler()
    
    for vol in volSet:
        index, fn = vol.getLocation()
        convert = True # by default, convert, which is the save way
        if index == em.NO_INDEX: # means single volumes
            volName = os.path.basename(fn)
            if fn.endswith('.em') or fn.endswith('.mrc'):
                convert = False # we can just create a link in this case
        else:
            volName = 'volume_%03d.mrc' % vol.getObjId()
        
        volFn = os.path.join(volDir, volName)
        if convert:
            ih.convert(vol, volFn)
        else:
            pwutils.createLink(fn, volFn)
        
        # Make the volumes names relative to the xml file
        # where the programs will be executed
        volRel = os.path.relpath(volFn, os.path.dirname(volXml))
        p = Particle()
        s = xcfScore()
        s.setValue(1.0)
        pytomInfo = getattr(vol, 'pytomInfo', None)
        if pytomInfo is None:
            p.setWedge(w)
        else:
            p.fromXML(pytomInfo.get()) # Get stored XML format from PyTom
        p.setFilename(volRel)
        p.setScore(s)
        pl.append(p)
        
    pl.toXMLFile(volXml)
    #pl.setWedgeAllParticles(w)
    sys.path = backupPath
    
    
def readSetOfVolumes(volsXml, volSet, **kwargs):
    """ Populate a Scipion set of volumes from a given
    xml file from PyTom.
    """
    # Add to the path the root to pytom
    backupPath = list(sys.path)
    addPyTomPaths()
    
    from pytom.basic.structures import Particle, ParticleList, Wedge
    from pytom.score.score import Score, PeakPrior
    from pytom.frm.FRMAlignment import FRMScore
    
    pl = ParticleList()
    pl.fromXMLFile(volsXml)
    
    xmlDir = os.path.dirname(volsXml)
    for particle in pl:
        volFn = os.path.join(xmlDir, particle.getFilename())
        vol = em.Volume()
        vol.setFileName(volFn)
        vol.pytomInfo = readPyTomInfo(particle)
        volSet.append(vol)
    
    sys.path = backupPath
    

def writeVolumesSqlite(volsXml, volsSqlite, **kwargs):
    """ Convert a volume list from PyTom xml format to Scipion sqlite file.
    If the volSqlite exists, it will be deleted.
    """
    pwutils.cleanPath(volsSqlite)
    volSet = em.SetOfVolumes(filename=volsSqlite)    
    readSetOfVolumes(volsXml, volSet, **kwargs)
    volSet.write()
    
    
def readPyTomInfo(particle):
    """ Extract all info from PyTom Particle object. """
    #NOTE: Right now the easy way is to store all information
    # in the PyTom native format that is xml, maybe later 
    # could be useful to standardize with the conventions.
    pytomInfo = pwobj.String()
    pytomInfo.set(str(particle))
    
    return pytomInfo
    
    
def readRotation(rotation):
    """ Convert a PyTom Rotation object to Scipion. """
    # Sample from XML:
    # <Rotation Paradigm="ZXZ" X="0.133725303909" Z1="359.345252413" Z2="0.260725774937"/>
    d = {'Paradigm': rotation.getParadigm(),
         'Z1': rotation.getZ1(),
         'X': rotation.getX(),
         'Z2': rotation.getZ2()}
    r = pwobj.String(value=json.dumps(d))
    
    return r


def readShift(shift):
    """ Convert from PyTom Shift object to Scipion. """
    #<Shift X="0.0293072589139" Y="0.0244929888746" Z="0.166532638949"/>
    #<PickPosition Origin="" X="0.0" Y="0.0" Z="0.0"/>
    d = {'X': shift.getX(),
         'Y': shift.getY(),
         'Z': shift.getZ()}
    s = pwobj.String(value=json.dumps(d))
    
    return s
    

def readPickPosition(pickPos):
    pass
    
def printLicense():
    """ Print the PyTom license file. """
    f = open(os.path.join(os.environ['PYTOM_HOME'], 'LICENSE.txt'))
    print "=================================================================="
    print "=                       PyTom LICENCE                            ="
    print "=================================================================="
    
    for line in f:
        sys.stdout.write(line)
    f.close()

