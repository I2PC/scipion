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

import collections
import os
from pyworkflow.em.convert import ImageHandler
import pyworkflow.utils as pwutils
import struct
import getpass

cootPdbTemplateFileName = "cootOut%04d.pdb"
cootScriptFileName = "cootScript.py"
ORIGIN = 0  # save coordinate origin in the mrc header field=Origin (Angstrom)
START = 1  # save coordinate origin in the mrc header field=start (pixel)


def getEnviron(ccp4First=True):
    environ = pwutils.Environ(os.environ)
    pos = pwutils.Environ.BEGIN if ccp4First else pwutils.Environ.END
    _ccp4_home = os.environ['CCP4_HOME']
    _ccp4_master, _dir = os.path.split(_ccp4_home)
    _username = getpass.getuser()
    environ.update({
            'PATH': os.path.join(_ccp4_home, 'bin'),
            'LD_LIBRARY_PATH': os.path.join(_ccp4_home, 'lib'),
            'CCP4_MASTER': _ccp4_master,
            'CCP4': _ccp4_home,
            'CCP4_SCR': os.path.join("/tmp", _username),
            'GFORTRAN_UNBUFFERED_PRECONNECTED':"Y",
            'CBIN': os.path.join(_ccp4_home, 'bin'),
            'CLIB': os.path.join(_ccp4_home, 'lib'),
            'CLIBD': os.path.join(_ccp4_home, 'lib/data'),
            'CETC': os.path.join(_ccp4_home, 'etc'),
            'CINCL': os.path.join(_ccp4_home, 'include'),
            'CHTML': os.path.join(_ccp4_home, 'html'),
            'CEXAM': os.path.join(_ccp4_home, 'examples'),
            'CCP4I_TOP': os.path.join(_ccp4_home, 'share/ccp4i'),
            'MMCIFDIC': os.path.join(_ccp4_home, 'lib', 'cif_mmdic.lib'),
            'CLIBD_MON': os.path.join(_ccp4_home, 'lib', 'data', 'monomers'),
            'CRANK': os.path.join(_ccp4_home, 'crank'),
            'CCP4_HELPDIR': os.path.join(_ccp4_home, 'help'),
            }, position=pos)
    return environ
#     """ DEFAULT CCP4 configuration file # This is CCP4 setup file for users of bash, zsh and similar shells.
# # The include/ccp4.setup-sh.in file is a template used to generate
# # the bin/ccp4.setup-sh file.  Source the latter to use CCP4.
# # You may add it to one of your shell startup files, such as ~/.bashrc:
# # . /full/path/to/ccp4/ccp4-version/bin/ccp4.setup-sh
#
# # To remove previously added CCP4 directories from the PATH, uncomment:
# #[ -n "$CCP4" ] && PATH="$(echo "$PATH" | sed "s,${CCP4}[^:]*:,,g")"
#
# ################### MOST-OFTEN CUSTOMIZED VARIABLES #####################
#
# # CCP4_MASTER is the location of the top-level directory containing ccp4-N.N.N.
# export CCP4_MASTER=/home/roberto
# export CCP4=$CCP4_MASTER/ccp4-6.5
#
# # CCP4_SCR: a per-user directory for run-time-generated scratch files.
# export CCP4_SCR=/tmp/`whoami`
#
# # CCP4I_TCLTK - directory containing tclsh, wish and bltwish executables
# #               for ccp4i. The CCP4 suite is distributed with Tcl/Tk 8.4.
# export CCP4I_TCLTK=$CCP4/bin
#
# # BALBES pipeline with large database
# test -d $CCP4_MASTER/BALBES && export BALBES_ROOT=$CCP4_MASTER/BALBES
#
# ### Optional - setting http proxy
# #export HTTP_PROXY=wwwblah.blah.ac.uk:xxxx/blah.blah
#
#
# ######################### ALL THE REST ###########################
#
# # This variable is set to ensure that the logfile output from programs
# # compiled with Gfortran is in the correct order.
# export GFORTRAN_UNBUFFERED_PRECONNECTED=Y
#
# # CBIN: location of the executables -- must be on your path (see below)
# export CBIN=$CCP4/bin
# # CLIB: location of (binary) library files such as libccp4.a and libccp4.so
# export CLIB=$CCP4/lib
# # CLIBD: platform-independent data files
# export CLIBD=$CCP4/lib/data
# # CETC: executable scripts (NOT configuration files)
# export CETC=$CCP4/etc
# # CINCL: headers and two *.def files for handling "logical names" in CCP4
# export CINCL=$CCP4/include
# # CHTML: html documentation
# export CHTML=$CCP4/html
# # CEXAM: examples and some tests
# export CEXAM=$CCP4/examples
# # source code directories
# #export CLIBS=$CCP4/lib/libccp4
# #export CPROG=$CCP4/src
# # CCP4I_TOP: the top directory of the interface
# export CCP4I_TOP=$CCP4/share/ccp4i
# # MMCIFDIC: platform-dependent (not in $CLIBD) data file for the ccif library
# export MMCIFDIC=$CLIB/ccp4/cif_mmdic.lib
# # CRANK: location of Crank automation suite within ccp4i
# export CRANK=$CCP4I_TOP/crank
# # CLIBD_MON: dictionary files for REFMAC5 (keep trailing /)
# export CLIBD_MON=$CCP4/lib/data/monomers/
# # CCP4_HELPDIR: location of the VMS-style help file used by (ip)mosflm
# export CCP4_HELPDIR=$CCP4/help/            # NB trailing /
#
# if test ! -d $CCP4; then
#     echo "WARNING: The directory $CCP4"
#     echo 'WARNING: (assigned to $CCP4) does not exist.'
#     echo "WARNING: The CCP4 programs will not run correctly."
# fi
#
# # check to see if $CCP4_SCR exists and if not try to make it
# test -d $CCP4_SCR || mkdir $CCP4_SCR
# test -d $CCP4_SCR || echo "Unable to assign CCP4_SCR. This will cause probs."
#
# # HARVESTHOME specifies location of harvesting files (defaults to $HOME)
# #export HARVESTHOME=$HOME
#
# # MOSFLM_WISH: if not set, ${CCP4I_TCLTK}/wish is used
# #export MOSFLM_WISH=${CCP4I_TCLTK}/wish
#
# # CCP4_OPEN is set to 'UNKNOWN' by default. When set to 'NEW' traditional
# # CCP4 programs treat over-writing of existing files as a fatal error.
# export CCP4_OPEN=UNKNOWN
#
# PATH=$CCP4/etc:$CCP4/bin:$CCP4/share/xia2/Applications:$PATH
#
# if [ -n "$MANPATH" ]; then
#   export MANPATH=$CCP4/share/man:$MANPATH
# fi
#
# # traditional aliases
# alias ccp4='pushd $CCP4>/dev/null'
# alias xtal='pushd $CCP4_MASTER>/dev/null'
# alias cbin='pushd $CBIN>/dev/null'
# alias cetc='pushd $CETC>/dev/null'
# #alias cprog='pushd $CPROG>/dev/null'
# alias cincl='pushd $CINCL>/dev/null'
# alias clib='pushd $CLIB>/dev/null'
# alias clibd='pushd $CLIBD>/dev/null'
# #alias clibs='pushd $CLIBS>/dev/null'
# alias cbin='pushd $CBIN>/dev/null'
# alias cexam='pushd $CEXAM>/dev/null'
# alias chtml='pushd $CHTML>/dev/null'
# """

def runCCP4Program(program, args="", extraEnvDict=None):
    """ Internal shortcut function to launch a CCP4 program. """
    env = getEnviron()
    # env.update(_envDict)
    if extraEnvDict is not None:
        env.update(extraEnvDict)
    pwutils.runJob(None, program, args, env=env)


def adaptFileToCCP4(inFileName, outFileName, scipionOriginShifts,
                    sampling=1.0, originField=START):
    """ Check input file format.
        if mrc, check if header and scipion database agree (regarding origin)
        if header and scipion object have the same origin creates a link to
        original file. Otherwise copy the file and fix the origin
    """
    # def compareFloatTuples(t1, t2, error=0.001):
    #    return abs(sum(tuple(x - y for x, y in zip(t1, t2)))) < error

    # scipion origin follows different conventions from ccp4
    scipionOriginShifts = tuple([-1. * x for x in scipionOriginShifts])
    x, y, z, ndim = ImageHandler().getDimensions(inFileName)
    # if inFileName.endswith('.mrc') or inFileName.endswith('.map')\
    #        or inFileName.endswith(':mrc'):
    #    ccp4header = Ccp4Header(inFileName, readHeader=True)
    #    ccp4Origin = ccp4header.getOrigin()
    #    if compareFloatTuples(ccp4Origin, scipionOriginShifts):
    #        if compareFloatTuples(ccp4header.getGridSampling(), (x, y, z)):
    #            if compareFloatTuples(ccp4header.getCellDimensions(),
    #                                  (x*sampling, y*sampling, z*sampling)):
    #                if inFileName.endswith(':mrc'):
    #                    pwutils.createLink(os.path.abspath(inFileName[:-4]),
    #                                    outFileName)
    #                return
    #
    # else:
    #    print "FILE %s does not end with mrc" % inFileName
    if z == 1:
        z = ndim

    ImageHandler().convert(inFileName, outFileName)
    ccp4header = Ccp4Header(outFileName, readHeader=True)
    ccp4header.setGridSampling(x, y, z)
    ccp4header.setCellDimensions(x * sampling, y * sampling, z * sampling)
    if originField == ORIGIN:
        ccp4header.setOrigin(scipionOriginShifts)
    else:
        ccp4header.setStartAngstrom(scipionOriginShifts, sampling)

    ccp4header.writeHeader()


def copyMRCHeader(inFileName, outFileName, scipionOriginShifts,
                  sampling, originField=START):
    x, y, z, ndim = ImageHandler().getDimensions(inFileName)
    ccp4header = Ccp4Header(outFileName, readHeader=True)
    ccp4header.setGridSampling(x, y, z)
    ccp4header.setCellDimensions(x * sampling, y * sampling, z * sampling)
    scipionOriginShifts = tuple([-1. * x for x in scipionOriginShifts])
    if originField == ORIGIN:
        ccp4header.setOrigin(scipionOriginShifts)
    else:
        ccp4header.setStartAngstrom(scipionOriginShifts, sampling)
    ccp4header.writeHeader()


def getProgram(progName):
    """ Return the program binary that will be used. """
    if 'CCP4_HOME' not in os.environ:
        return None
    return os.path.join(os.environ['CCP4_HOME'], 'bin',
                        os.path.basename(progName))


class Ccp4Header():
    """
    In spite of the name this is the MRC2014 format no the CCP4.
    In an MRC EM file the origin is typically in header fields called
    xorigin, yorigin, zorigin which are specific to EM data (MRC 2000 format)
    and are not part of the nearly identical CCP4 format used for x-ray maps
    that Coot is expecting.
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
      17      MAPC            Which axis corresponds to Cols. (1,2,3 for X,Y,Z)
      18      MAPR            Which axis corresponds to Rows. (1,2,3 for X,Y,Z)
      19      MAPS            Which axis corresponds to Sects.
                              (1,2,3 for X,Y,Z)
      20      AMIN            Minimum density value
      21      AMAX            Maximum density value
      22      AMEAN           Mean    density value    (Average)
      23      ISPG            Space group number
      24      NSYMBT          Number of bytes used for storing symmetry
                              operators
      25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
      26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc)
                              if LSKFLG .ne. 0.
      35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
                              Skew transformation is from standard orthogonal
                              coordinate frame (as used for atoms)
                              to orthogonal map frame, as
                              Xo(map) = S * (Xo(atoms) - t)
      38      future use      (some of these are used by the MSUBSX routines
       .          "           in MAPBRICK, MAPCONT and FRODO)
       .          "           (all set to zero by default)
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

    def __init__(self, fileName, readHeader=False):
        self._name = fileName.replace(':mrc', '')  # remove mrc ending
        self._header = collections.OrderedDict()
        self.chain = "< 3i i 3i 3i 3f 144s 3f"

        if readHeader:
            self.readHeader()

    def setOrigin(self, originTransformShift):
        # TODO: should we use originX,Y,Z and set this to 0
        self._header['originX'] = originTransformShift[0]
        self._header['originY'] = originTransformShift[1]
        self._header['originZ'] = originTransformShift[2]
        self._header['NCSTART'] = 0.
        self._header['NRSTART'] = 0.
        self._header['NSSTART'] = 0.

    def getOrigin(self):
        return self._header['originX'], \
               self._header['originY'], \
               self._header['originZ']

    def setSampling(self, sampling):
        self._header['Xlength'] = self._header['NX'] * sampling
        self._header['Ylength'] = self._header['NY'] * sampling
        self._header['Zlength'] = self._header['NZ'] * sampling

    def setStartPixel(self, originTransformShift):  # PIXEL
        """input pixels"""
        self._header['originX'] = 0.  # originTransformShift[0]
        self._header['originY'] = 0.  # originTransformShift[1]
        self._header['originZ'] = 0.  # originTransformShift[2]
        self._header['NCSTART'] = originTransformShift[0]
        self._header['NRSTART'] = originTransformShift[1]
        self._header['NSSTART'] = originTransformShift[2]

    def setStartAngstrom(self,  originTransformShift, sampling):  # Angstrom
        """input Angstrom"""
        self.setStartPixel(tuple([int(round(x/sampling)) for x in
                           originTransformShift]))

    def getStartPixel(self):  # PIXEL
        """returns pixels"""
        return self._header['NCSTART'],\
               self._header['NRSTART'], \
               self._header['NSSTART']

    def getStartAngstrom(self, sampling):  # Angstrom
        """returns Angstrom"""
        return tuple([x*sampling for x in self.getStartPixel()])

    def getDims(self):
        return self._header['NC'],\
               self._header['NR'],\
               self._header['NS']

    def setGridSampling(self, x, y, z):
        self._header['NX'] = x
        self._header['NY'] = y
        self._header['NZ'] = z

    def getGridSampling(self):
        return self._header['NX'],\
               self._header['NY'],\
               self._header['NZ']

    def setCellDimensions(self, x, y, z):
        self._header['Xlength'] = x
        self._header['Ylength'] = y
        self._header['Zlength'] = z

    def getCellDimensions(self):
        return self._header['Xlength'],\
               self._header['Ylength'],\
               self._header['Zlength']

    def readHeader(self):
        # check file exists

        # read header
        f = open(self._name, 'rb')
        s = f.read(52*4)  # read header from word 0 to 51
        f.close()
        a = struct.unpack(self.chain, s)

        # fill dicionary
        # self._header['empty']=""
        self._header['NC'] = a[0]
        self._header['NR'] = a[1]
        self._header['NS'] = a[2]
        self._header['Mode'] = a[3]
        self._header['NCSTART'] = a[4]
        self._header['NRSTART'] = a[5]
        self._header['NSSTART'] = a[6]
        self._header['NX'] = a[7]
        self._header['NY'] = a[8]
        self._header['NZ'] = a[9]
        self._header['Xlength'] = a[10]
        self._header['Ylength'] = a[11]
        self._header['Zlength'] = a[12]
        self._header['dummy1'] = a[13] + "\n"  # "< 3i i 3i 3i 3f 36s 3f"
        self._header['originX'] = a[14]
        self._header['originY'] = a[15]
        self._header['originZ'] = a[16]

    def getHeader(self):
        return self._header

    def setHeader(self, newHeader):
        self._header = newHeader

    def writeHeader(self):
        ss = struct.Struct(self.chain)
        t = tuple(self._header.values())
        packed_data = ss.pack(*t)
        f = open(self._name, 'r+')
        f.write(packed_data)
        f.close()

    def __str__(self):
        s = ""
        for k, v in self._header.iteritems():
            s += "%s: %s\n" % (str(k), str(v))
        return s

    def computeSampling(self):
        return self._header['Zlength'] / self._header['NZ']
