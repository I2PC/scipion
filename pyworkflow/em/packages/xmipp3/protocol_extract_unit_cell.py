# **************************************************************************
# *
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
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import StringParam, PointerParam, FloatParam, EnumParam
from pyworkflow.em.constants import SYM_I222

class XmippProtExtractUnit(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'extract unit cell'
    _program = "" 
    _version = VERSION_1_2

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # def _createFileNamesTemplates(self):
    #     myDict = {
    #               'volume' : 'final_volume.mrc',
    #               'fsc' : 'final_fsc.xml'
    #               }
    #     self._updateFilenamesDict(myDict)

        #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, label="Input Volume", important=True,
                      pointerClass='Volume',
                      help='This volume will be cropped')
        form.addParam('symmetryGroup', EnumParam, choices=["I2 (I222)"],
                      default=SYM_I222,
                      label="Symmetry",
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry"
                           " for a description of the symmetry groups format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                      )
        #WIZARD
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
        self._createFileNamesTemplates()
        self._insertFunctionStep('exportVolumeStep')
        self._insertFunctionStep('exportFSCStep')

    #--------------------------- STEPS functions --------------------------------------------

    def exportVolumeStep(self):

        ih = ImageHandler()
        ih.convert(self.exportVolume.get().getLocation(), self.getFnPath())

    def exportFSCStep(self):

        x,y = self.exportFSC.get().getData()
        fo = open(self.getFnPath("fsc"), "w")
        fo.write('<fsc title="FSC(%s)" xaxis="Resolution(A-1)" '
                 'yaxis="Correlation Coefficient">\n' % self._getFileName('volume'))
        for i in range(len(x)):
            fo.write("<coordinate>\n")
            fo.write("<x>%f</x>\n"%x[i])
            fo.write("<y>%f</y>\n" % y[i])
            fo.write("</coordinate>\n")

        fo.write("</fsc>\n")
        fo.close()

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

#    def getFnPath(self, label='volume'):
#        return os.path.join(self.filesPath.get(), self._getFileName(label))
