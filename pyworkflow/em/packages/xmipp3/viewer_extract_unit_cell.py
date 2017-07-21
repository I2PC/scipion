# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es), May 2013
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

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer
from pyworkflow.em.packages.chimera.convert import runChimeraProgram, getProgram, symMapperScipionchimera
from protocol_extract_unit_cell import XmippProtExtractUnit
from pyworkflow.em.constants import SYM_I222
import os

class viewerXmippProtExtractUnit(XmippViewer):
    """ Visualize the output of protocol volume strain """
    _label = 'viewer extract unit cell'
    _targets = [XmippProtExtractUnit]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        XmippViewer.__init__(self, **args)

#"""ROB: I know that there is a nice chimera interface
# but it does not work in this case since I am interested
# in reading the MRC header. So I will use chimera as
# an external program"""

    def _visualize(self, obj, **args):
        #save temporal file
        sampling = self.protocol.inputVolumes.get().getSamplingRate()
        tmpFileName = self.protocol._getTmpPath("chimera.cmd")
        f = open(tmpFileName,"w")
        if self.protocol.symmetryGroup >= SYM_I222:
            cMap=['red', 'yellow', 'green', 'cyan', 'blue']
            d={}
            d['outerRadius']=self.protocol.outerRadius.get() * sampling
            d['innerRadius']=self.protocol.innerRadius.get() * sampling
            d['symmetry']=symMapperScipionchimera[self.protocol.symmetryGroup.get()]
            f.write("shape icosahedron mesh true radius %(outerRadius)d orientation %(symmetry)s\n"%d)
            step = (d['outerRadius'] - d['innerRadius']) / float(len(cMap)-1)
            f.write("scolor #1  geom radial center 0,0,0 cmap ")
            counter = 0
            s=""
            for color in cMap:
                s += "%d,%s:"%(d['innerRadius'] + counter * step, cMap[counter])
                counter += 1
            f.write(s[:-1]+'\n')

        volSize0 = self.protocol.inputVolumes.get().getDim()[0]
        origin0 = - sampling * volSize0 / 2.0
        origin1 =  self.protocol.outputVolume.getOrigin().getShifts()
        f.write("volume #0 style mesh    level .001 origin %f\n"%origin0)
        x = origin1[0] * sampling
        y = origin1[1] * sampling
        z = origin1[2] * sampling

        f.write("volume #1 style surface level .001 origin %f,%f,%f\n"%(x,y,z))
        f.write("center\n"  )
        f.close()

        inputVol  = self.protocol.inputVolumes.get().getFileName().replace(':mrc', '')
        outputVol = self.protocol.outputVolume.getFileName()
        args = " "
        if os.path.exists(inputVol):
            args += inputVol + " "
        if os.path.exists(outputVol):
            args += outputVol + " "
        if True and os.path.exists(tmpFileName):
            args += tmpFileName + " "
        runChimeraProgram(getProgram(), args)
