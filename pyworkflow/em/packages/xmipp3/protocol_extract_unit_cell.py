# **************************************************************************
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
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


from pyworkflow import VERSION_1_2
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import  PointerParam, FloatParam, EnumParam
from pyworkflow.em.constants import SYM_I222r
from pyworkflow.em.packages.xmipp3 import XMIPP_SYM_NAME
from pyworkflow.em.constants import SCIPION_SYM_NAME
from pyworkflow.em import Volume

class XmippProtExtractUnit(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'extract unit cell'
    _program = "" 
    _version = VERSION_1_2

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, label="Input Volume", important=True,
                      pointerClass='Volume',
                      help='This volume will be cropped')
        form.addParam('symmetryGroup', EnumParam, choices=[XMIPP_SYM_NAME[SYM_I222r] +
                                                           " (" + SCIPION_SYM_NAME[SYM_I222r] + ")"],
                      default=SYM_I222r,
                      label="Symmetry",
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry"
                           " for a description of the symmetry groups format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                      )
#form.addParam('symmetry', TextParam, default='C1',
#                      label='Point group symmetry:',
#                      condition='not doContinue',
#                      help='Parameter *ASYM* in FREALIGN\n\n'
#                           'Specify the symmetry.Choices are: Cn,Dn,T,O,I,I1,I2,N or H (can be zero)\n'
#                           'n  = rotational symmetry required in pointgroup C(n) or D(n)\n'
#                           'N  = number of symmetry matrices to read in.\n'
#                           'T  = tetrahedral pointgroup 23\n'
#                           'O  = octahedral pointgroup 432\n'
#                           'I  = icosahedral 532 symmetry in setting 1 (5fold is on X)\n'
#                           'I1 = also in setting 1 (X) - as used by Imagic\n'
#                           'I2 = in setting 2 (Y) - as used by Crowther et. al\n'
#~/scipion/pyworkflow/em/packages/grigoriefflab/protocol_frealign_base.py

        form.addParam('innerRadius', FloatParam, default=-1,
                      label="Inner Radius (px)", help="inner Mask radius, if -1, the radius will be 0")
        form.addParam('outerRadius', FloatParam, default=-1,
                      label="Outer Radius (px)", help="outer Mask radius, if -1, the radius will be volume_size/2")
        form.addParam('expandFactor', FloatParam, default=0.,
                      label="Expand Factor", help="Increment cropped region by this factor")
        form.addParam('offset', FloatParam, default=0.,
                      label="offset", help="to be defined")

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('extractUnit')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------

    def extractUnit(self):
        #        samplingRate = protocol._getSetSampling()
        args = "-i %s -o %s" % (self.inputVolumes.get().getFileName(), self._getOutputVol())
        args += " --unitcell %s "% XMIPP_SYM_NAME[self.symmetryGroup.get()]
        args += " %f "% self.innerRadius.get()
        args += " %f "% self.outerRadius.get()
        args += " %f "% self.expandFactor.get()
        args += " %f "% self.offset.get()
        print "args", args
        self.runJob("xmipp_transform_window", args)

    def createOutputStep(self):
        vol = Volume()
        vol.setLocation(self._getOutputVol())
        vol.setSamplingRate(self.inputVolumes.get().getSamplingRate())
        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputVolumes, self.outputVolume)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        message = []
        return message

    def _summary(self):
        #message = "Data Available at : *%s*"% self.filesPath.get()
        message=""
        return [message]

    def _methods(self):
        return []

#--------------------------- UTILS functions ---------------------------------------------------

def _getOutputVol(self):
    return self._getExtraPath("output_volume.mrc")