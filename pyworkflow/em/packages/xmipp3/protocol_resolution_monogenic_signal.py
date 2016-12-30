# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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



from pyworkflow.protocol.params import (PointerParam, StringParam, BooleanParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtRefine3D
from convert import readSetOfVolumes
from pyworkflow.em import ImageHandler
from pyworkflow.utils import getExt, removeExt
from os.path import abspath
import numpy as np

OUTPUT_RESOLUTION_FILE = 'mgresolution.vol'
FN_FILTERED_MAP = 'filteredMap.vol'
OUTPUT_RESOLUTION_FILE_CHIMERA = 'MG_Chimera_resolution.vol'
FN_MEAN_VOL = 'mean_volume.vol'


class XmippProtMonoRes(ProtRefine3D):
    """    
    Given a map the protocol assigns local resolutions to each pixel of the map.
    """
    _label = 'monogenic resolution'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('halfVolumes', BooleanParam, default=False,
                      label="Would you like to use half volumes?",
                      help='The noise estimation for determining the local resolution '
                           'is performed via half volumes.')

        form.addParam('inputVolumes', PointerParam, pointerClass='Volume',
                      label="Input Volume",
                      condition = 'not halfVolumes',
                      help='Select a volume for determining its local resolution.')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input Volume",
                      condition = 'halfVolumes',
                      help='Select a volume for determining its local resolution.')

        form.addParam('inputVolume2', PointerParam, pointerClass='Volume',
                      label="Second Half Volume",
                      condition='halfVolumes',
                      help='Select a second volume for determining a local resolution.')

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask', 
                      condition='(halfVolumes) or (not halfVolumes)',
                      label="Binary Mask",
                      help='The mask determines which points are specimen and which ones not')

        form.addParam('symmetry', StringParam, default='c1',
                      label="Symmetry",
                      help='Symmetry group. By default = c1')

        line = form.addLine('Resolution Range (A)',
                            help="If the user knows the range of resolutions or only a"
                                 " range of frequency needs to be analysed", expertLevel=LEVEL_ADVANCED)
        form.addParam('isPremasked', BooleanParam, default=False,
                      label="Is the original premasked?",
                      help='Sometimes the original volume is masked inside a spherical mask. In this case'
                      'please select yes')
        
        form.addParam('volumeRadius', FloatParam, default=-1,
                      label="Spherical mask radius",
                      condition = 'isPremasked and not halfVolumes', 
                      help='When the origianl volume is originally premasked, the noise estimation ought'
                      'to be performed inside that premask, and out of the provieded mask asked in the previus'
                      'box. The radius value, determines the radius of the spherical premask. By default'
                      'radius = -1 use the half of the volume size as radius')
        
        form.addParam('volumeRadiusHalf', FloatParam, default=-1,
                      label="Spherical mask radius",
                      condition = 'halfVolumes and isPremasked',
                      help='When the origianl volume is originally premasked, the noise estimation ought'
                      'to be performed inside that premask, and out of the provieded mask asked in the previus'
                      'box. The radius value, determines the radius of the spherical premask. By default'
                      'radius = -1 use the half of the volume size as radius')

        line.addParam('minRes', FloatParam, default=1, label='High')
        line.addParam('maxRes', FloatParam, default=50, label='Low')

        form.addParam('filterInput', BooleanParam, default=False, 
                      label="Filter input volume with local resolution?",
                      help='The input map is locally filtered at the local resolution map.')

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        
        self.micsFn = self._getPath()
        if (not self.halfVolumes):
            self.vol0Fn = self.inputVolumes.get().getFileName()
            self.maskFn = self.Mask.get().getFileName()

        if self.halfVolumes.get() is True:
            self.vol1Fn = self.inputVolume.get().getFileName()
            self.vol2Fn = self.inputVolume2.get().getFileName()
            self.maskFn = self.Mask.get().getFileName()

            # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep', )

        MS = self._insertFunctionStep('resolutionMonogenicSignalStep', prerequisites=[convertId])

        self._insertFunctionStep('createOutputStep', prerequisites=[MS])

        self._insertFunctionStep("createChimeraScript")

        self._insertFunctionStep("createHistrogram")

    def convertInputStep(self):
        """ Read the input volume.
        """
        if (self.halfVolumes.get() is False):
            extVol0 = getExt(self.vol0Fn)
            if (extVol0 == '.mrc') or (extVol0 == '.map'):
                self.vol0Fn = self.vol0Fn + ':mrc'

        if self.halfVolumes.get() is True:
            extVol1 = getExt(self.vol1Fn)
            extVol2 = getExt(self.vol2Fn)
            if (extVol1 == '.mrc') or (extVol1 == '.map'):
                self.vol1Fn = self.vol1Fn + ':mrc'
            if (extVol2 == '.mrc') or (extVol2 == '.map'):
                self.vol2Fn = self.vol2Fn + ':mrc'

        extMask = getExt(self.maskFn)
        if (extMask == '.mrc') or (extMask == '.map'):
            self.maskFn = self.maskFn + ':mrc'

    def resolutionMonogenicSignalStep(self):

        if (self.halfVolumes):
            if (self.isPremasked):
                if (self.volumeRadiusHalf.get() is -1):
                    xdim, _ydim, _zdim = self.inputVolume.get().getDim()
                    xdim = xdim*0.5
                else:
                    xdim = self.volumeRadiusHalf.get()
            else:
                xdim, _ydim, _zdim = self.inputVolume.get().getDim()
                xdim = xdim*0.5
        else:
            if (self.isPremasked):
                if (self.volumeRadius.get() is -1):
                    xdim, _ydim, _zdim = self.inputVolumes.get().getDim()
                    xdim = xdim*0.5
                else:
                    xdim = self.volumeRadius.get()
            else:
                xdim, _ydim, _zdim = self.inputVolumes.get().getDim()
                xdim = xdim*0.5
                
        if self.halfVolumes.get() is False:
            params = ' --vol %s' % self.vol0Fn
            params += ' --mask %s' % self.maskFn
        else:
            params = ' --vol %s' % self.vol1Fn
            params += ' --vol2 %s' % self.vol2Fn
            params += ' --meanVol %s' % self._getExtraPath(FN_MEAN_VOL)
            params += ' --mask %s' % self.maskFn

        params += ' -o %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        if (self.halfVolumes):
            params += ' --sampling_rate %f' % self.inputVolume.get().getSamplingRate()
        else:
            params += ' --sampling_rate %f' % self.inputVolumes.get().getSamplingRate()
        params += ' --number_frequencies %f' % 50
        params += ' --minRes %f' % self.minRes.get()
        params += ' --maxRes %f' % self.maxRes.get()
        params += ' --volumeRadius %f' % xdim
        params += ' --chimera_volume %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE_CHIMERA)
        params += ' --linear '
        params += ' --sym %s' % self.symmetry.get()
        params += ' --trimmed %f' % 98  #This parameter only considers resolution values in percentile 98
        if self.filterInput.get():
            params += ' --filtered_volume %s' % self._getExtraPath(FN_FILTERED_MAP)
        else:
            params += ' --filtered_volume %s' % ''



        self.runJob('xmipp_resolution_monogenic_signal', params)

    def createChimeraScript(self):
        fnRoot = "extra/"
        scriptFile = self._getPath('Chimera_resolution.cmd')
        fhCmd = open(scriptFile, 'w')
        imageFile = self._getExtraPath(OUTPUT_RESOLUTION_FILE_CHIMERA)
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        min_Res = round(np.amin(imgData)*100)/100
        max_Res = round(np.amax(imgData)*100)/100
        inter = round((max_Res - min_Res)*20)/100
        print inter
        
        if self.halfVolumes.get() is True:
            #fhCmd.write("open %s\n" % (fnRoot+FN_MEAN_VOL)) #Perhaps to check the use of mean volume is useful
            fnbase = removeExt(self.inputVolume.get().getFileName())
            ext = getExt(self.inputVolume.get().getFileName())
            fninput = abspath(fnbase + ext[0:4])
            fhCmd.write("open %s\n" % fninput)
        else:
            fnbase = removeExt(self.inputVolumes.get().getFileName())
            ext = getExt(self.inputVolumes.get().getFileName())
            fninput = abspath(fnbase + ext[0:4])
            fhCmd.write("open %s\n" % fninput)
        fhCmd.write("open %s\n" % (fnRoot + OUTPUT_RESOLUTION_FILE_CHIMERA))
        if self.halfVolumes.get() is True:
            smprt = self.inputVolume.get().getSamplingRate()
        else:
            smprt = self.inputVolumes.get().getSamplingRate()
        fhCmd.write("volume #1 voxelSize %s\n" % (str(smprt)))
        fhCmd.write("vol #1 hide\n")
        fhCmd.write("scolor #0 volume #1 cmap rainbow reverseColors True\n")
        fhCmd.write("colorkey 0,0 0.2,0.6 %s blue %s cyan %s green %s yellow %s red\n" %(min_Res, 
                                                                                         min_Res+inter, 
                                                                                         min_Res+2*inter,
                                                                                         min_Res+3*inter,
                                                                                         max_Res))       
        fhCmd.close()


    def createHistrogram(self):

        params = ' -i %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        params += ' --mask binary_file %s' % self.maskFn
        params += ' --steps %f' % 30
        params += ' --range %f %f' % (self.minRes.get(), self.maxRes.get())
        params += ' -o %s' % self._getExtraPath('hist.xmd')

        self.runJob('xmipp_image_histogram', params)

    def createOutputStep(self):
        volume_path = self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        self.volumesSet = self._createSetOfVolumes('resolutionVol')
        if (self.halfVolumes):
            self.volumesSet.setSamplingRate(self.inputVolume.get().getSamplingRate())
        else:
            self.volumesSet.setSamplingRate(self.inputVolumes.get().getSamplingRate())
        readSetOfVolumes(volume_path, self.volumesSet)
        self._defineOutputs(outputVolume=self.volumesSet)
        if (self.halfVolumes):
            self._defineSourceRelation(self.inputVolume, self.volumesSet)
        else:
            self._defineSourceRelation(self.inputVolumes, self.volumesSet)
        

        if self.filterInput.get():
            print 'Saving filtered map'
            volume_filtered_path = self._getExtraPath(FN_FILTERED_MAP)
            self.volumesSet2 = self._createSetOfVolumes('filteredVol')
            if (self.halfVolumes):
                self.volumesSet2.setSamplingRate(self.inputVolume.get().getSamplingRate())
            else:
                self.volumesSet2.setSamplingRate(self.inputVolumes.get().getSamplingRate())
            readSetOfVolumes(volume_filtered_path, self.volumesSet2)
            self._defineOutputs(outputVolume_Filtered=self.volumesSet2)
            if (self.halfVolumes):
                self._defineSourceRelation(self.inputVolume, self.volumesSet2)
            else:
                self._defineSourceRelation(self.inputVolumes, self.volumesSet2)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):

        validateMsgs = []
        if self.halfVolumes.get() is True:
            if (not self.inputVolume.get().hasValue()):
                validateMsgs.append('Please provide input volume.')
        else:
            if (not self.inputVolumes.get().hasValue()):
                validateMsgs.append('Please provide input volume.')
        return validateMsgs

    def _methods(self):
        messages = []
        if hasattr(self, 'outputVolume'):
            messages.append(
                'The local resolution is performed [Publication: Not yet]')
        return messages

    def _citations(self):
        return ['Not yet']
