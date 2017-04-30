# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
import pyworkflow.em as em
import pyworkflow.utils as pwutils
import struct


#
cootPdbTemplateFileName="scipionOut%04d.pdb"
cootScriptFileName = "cootScript.py"


def getEnviron(ccp4First=True):
    environ = pwutils.Environ(os.environ)
    pos = pwutils.Environ.BEGIN if ccp4First else pwutils.Environ.END
    environ.update({
            'PATH': os.path.join(os.environ['CCP4_HOME'], 'bin'),
            'LD_LIBRARY_PATH': os.path.join(os.environ['CCP4_HOME'], 'lib'),
            }, position=pos)
    return environ

def runCCP4Program(program, args="", extraEnvDict=None):
    """ Internal shortcut function to launch a CCP4 program. """
    env=getEnviron()
    #env.update(_envDict)
    print "extraEnvDict", extraEnvDict
    if extraEnvDict is not None:
        env.update(extraEnvDict)
    pwutils.runJob(None, program, args, env=env)

def adapBinFileToCCP4(inFileName,outFileName):
    if inFileName.endswith('.mrc'):
        pwutils.createLink(inFileName, outFileName)
    else:
        em.ImageHandler().convert(inFileName, outFileName)

def getProgram(progName):
    """ Return the program binary that will be used. """
    if (not 'CCP4_HOME' in os.environ):
        return None
    return os.path.join(os.environ['CCP4_HOME'], 'bin',
                           os.path.basename(progName))

class ccp4Header():
    """
        The header is organised as 56 words followed by space for ten 80
    character text labels as follows::
       1      NC              # of Columns    (fastest changing in map)
       2      NR              # of Rows
       3      NS              # of Sections   (slowest changing in map)
       4      MODE            Data type
                                0 = envelope stored as signed bytes (from
                                    -128 lowest to 127 highest)
                                1 = Image     stored as Integer*2
                                2 = Image     stored as Reals
                                3 = Transform stored as Complex Integer*2
                                4 = Transform stored as Complex Reals
                                5 == 0
                                Note: Mode 2 is the normal mode used in
                                      the CCP4 programs. Other modes than 2 and 0
                                      may NOT WORK
       5      NCSTART         Number of first COLUMN  in map
       6      NRSTART         Number of first ROW     in map
       7      NSSTART         Number of first SECTION in map
       8      NX              Number of intervals along X
       9      NY              Number of intervals along Y
      10      NZ              Number of intervals along Z
      11      X length        Cell Dimensions (Angstroms)
      12      Y length                     "
      13      Z length                     "
      14      Alpha           Cell Angles     (Degrees)
      15      Beta                         "
      16      Gamma                        "
      17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
      18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
      19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
      20      AMIN            Minimum density value
      21      AMAX            Maximum density value
      22      AMEAN           Mean    density value    (Average)
      23      ISPG            Space group number
      24      NSYMBT          Number of bytes used for storing symmetry operators
      25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
      26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
                              LSKFLG .ne. 0.
      35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
                              Skew transformation is from standard orthogonal
                              coordinate frame (as used for atoms) to orthogonal
                              map frame, as
                                      Xo(map) = S * (Xo(atoms) - t)
      38      future use       (some of these are used by the MSUBSX routines
       .          "              in MAPBRICK, MAPCONT and FRODO)
       .          "   (all set to zero by default)
       .          "
      52          "
      53    MAP             Character string 'MAP ' to identify file type
      54    MACHST          Machine stamp indicating the machine type
                              which wrote file
      55      ARMS            Rms deviation of map from mean density
      56      NLABL           Number of labels being used
      57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)
    Symmetry records follow - if any - stored as text as in International
    Tables, operators separated by ``*`` and grouped into 'lines' of 80
    characters (i.e. symmetry operators do not cross the ends of the
    80-character 'lines' and the 'lines' do not terminate in a ``*``).
    """

    def __init__(self, fileName, readHeader= False):
        self._name   = fileName.replace(':mrc', '')# remove mrc ending
        self._header = None
        self.chain = "< 3i i 3i 3i 3f 36i 3f"

        if readHeader:
            self.readHeader()

    def readHeader(self):
        #check file exists

        #read header
        f = open(self.fileName,'rb')
        s = f.read(13*4)#read header up to angles incluse word 6
        f.close()
        a = struct.unpack(self.chain, s)

        #fill dicionary
        self.header['empty']=""
        self.header['NC'] = a[0]
        self.header['NR'] = a[1]
        self.header['NS'] = a[2]
        self.header['Mode'] = a[3]
        self.header['NCSTART'] = a[4]
        self.header['NRSTART'] = a[5]
        self.header['NSSTART'] = a[6]
        self.header['NX'] = a[7]
        self.header['NY'] = a[8]
        self.header['NZ'] = a[9]
        self.header['Xlength'] = a[10]
        self.header['Ylength'] = a[11]
        self.header['Zlength'] = a[12]
        self.header['originX'] = a[49]
        self.header['originY'] = a[50]
        self.header['originZ'] = a[51]
        self.header['samplingRateZ'] = a[12]/a[9]

    def getHeader(self):
        return self._header

    def setHeader(self, newHeader):
        self._header = newHeader

    def writeHeader(self):
        pass
        ss=struct.Struct(self.chain)
        t=tuple(self._header.values())
        packed_data = ss.pack(*t)
        f = open(self._name, 'r+')
        f.write(packed_data)
        f.close()
