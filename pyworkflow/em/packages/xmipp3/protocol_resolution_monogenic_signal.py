# -*- coding: utf-8 -*-
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
from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import (PointerParam, StringParam, BooleanParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from convert import readSetOfVolumes
from pyworkflow.object import Float
from pyworkflow.em import ImageHandler
from pyworkflow.utils import getExt
import numpy as np
import pyworkflow.em.metadata as md


MONORES_METHOD_URL = 'http://github.com/I2PC/scipion/wiki/XmippProtMonoRes'

OUTPUT_RESOLUTION_FILE = 'mgresolution.vol'
FN_FILTERED_MAP = 'filteredMap.vol'
OUTPUT_RESOLUTION_FILE_CHIMERA = 'MG_Chimera_resolution.vol'
OUTPUT_MASK_FILE = 'output_Mask.vol'
FN_MEAN_VOL = 'mean_volume.vol'
METADATA_MASK_FILE = 'mask_data.xmd'


class XmippProtMonoRes(ProtAnalysis3D):
    """    
    Given a map the protocol assigns local resolutions to each voxel of the map.
    """
    _label = 'local MonoRes'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.min_res_init = Float() 
        self.max_res_init = Float() 
        
    
    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('halfVolumes', BooleanParam, default=False,
                      label="Would you like to use half volumes?",
                      help='The noise estimation for determining the local resolution '
                           'is performed via half volumes.')

        form.addParam('inputVolumes', PointerParam, pointerClass='Volume',
                      label="Input Volume", important=True,
                      condition = 'not halfVolumes',
                      help='Select a volume for determining its local resolution.')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Volume Half 1", important=True,
                      condition = 'halfVolumes', 
                      help='Select a volume for determining its local resolution.')

        form.addParam('inputVolume2', PointerParam, pointerClass='Volume',
                      label="Volume Half 2", important=True,
                      condition='halfVolumes',
                      help='Select a second volume for determining a local resolution.')

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask', 
                      condition='(halfVolumes) or (not halfVolumes)',
                      label="Binary Mask", important=True,
                      help='The mask determines which points are specimen and which ones not')

        group = form.addGroup('Extra parameters')
        group.addParam('symmetry', StringParam, default='c1',
                      label="Symmetry",
                      help='Symmetry group. By default = c1.'
                      'See [[http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry][Symmetry]]'
                      'for a description of the symmetry groups format, If no symmetry is present, give c1.')

        line = group.addLine('Resolution Range (A)',
                            help="If the user knows the range of resolutions or only a"
                                 " range of frequency needs to be analysed")
        
        group.addParam('significance', FloatParam, default=0.95, expertLevel=LEVEL_ADVANCED,
                      label="Significance",
                      help='Relution is computed using hipothesis tests, this value determines'
                      'the significance of that test')
        
        group.addParam('isPremasked', BooleanParam, default=False,
                      label="Is the original premasked?",
                      help='Sometimes the original volume is masked inside a spherical mask. In this case'
                      'please select yes')
        
        form.addParam('noiseonlyinhalves', BooleanParam, expertLevel=LEVEL_ADVANCED,
                      default=True,
                      label="Use noise outside the mask?", 
                      condition = 'halfVolumes',
                      help='Select yes if the volume present noise outside the mask.'
                      ' Otherwise, select No.')
        
        group.addParam('volumeRadius', FloatParam, default=-1,
                      label="Spherical mask radius (px)",
                      condition = 'isPremasked and not halfVolumes', 
                      help='When the original volume is originally premasked, the noise estimation ought'
                      'to be performed inside that premask, and out of the provieded mask asked in the previus'
                      'box. The radius value, determines the radius of the spherical premask. By default'
                      'radius = -1 use the half of the volume size as radius')
        
        group.addParam('volumeRadiusHalf', FloatParam, default=-1,
                      label="Spherical mask radius (px)",
                      condition = 'halfVolumes and isPremasked',
                      help='When the origianl volume is originally premasked,'
                      'the noise estimation ought to be performed inside that'
                      'premask, and out of the provieded mask asked in the previus'
                      'box. The radius value, determines the radius in pixels of '
                      'the spherical premask. By default radius = -1 use the half'
                      'of the volume size as radius')

        line.addParam('minRes', FloatParam, default=1, label='High')
        line.addParam('maxRes', FloatParam, default=30, label='Low')
        line.addParam('stepSize', FloatParam, allowsNull=True,
                      expertLevel=LEVEL_ADVANCED, label='Step')

        group.addParam('filterInput', BooleanParam, default=False, 
                      label="Filter input volume with local resolution?",
                      help='The input map is locally filtered at the local resolution map.')

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        
        self.micsFn = self._getPath()

        if self.halfVolumes:
            self.vol1Fn = self.inputVolume.get().getFileName()
            self.vol2Fn = self.inputVolume2.get().getFileName()
            self.maskFn = self.Mask.get().getFileName()

            self.inputVolumes.set(None)

        else:
            self.vol0Fn = self.inputVolumes.get().getFileName()
            self.maskFn = self.Mask.get().getFileName()
            self.inputVolume.set(None)
            self.inputVolume2.set(None)

            # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep', )

        MS = self._insertFunctionStep('resolutionMonogenicSignalStep',
                                      prerequisites=[convertId])

        self._insertFunctionStep('createOutputStep', prerequisites=[MS])

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

        # Number of frequencies
        if self.stepSize.hasValue():
            Nfreqs = round((self.maxRes.get() - self.minRes.get())/self.stepSize.get())
        else:
            Nfreqs = 50
  
        if self.halfVolumes:
            if self.isPremasked:
                if self.volumeRadiusHalf == -1:
                    xdim, _ydim, _zdim = self.inputVolume.get().getDim()
                    xdim = xdim*0.5
                else:
                    xdim = self.volumeRadiusHalf.get()
            else:
                xdim, _ydim, _zdim = self.inputVolume.get().getDim()
                xdim = xdim*0.5
        else:
            if self.isPremasked:
                if self.volumeRadius == -1:
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
        params += ' --mask_out %s' % self._getExtraPath(OUTPUT_MASK_FILE)
        params += ' -o %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        if (self.halfVolumes):
            params += ' --sampling_rate %f' % self.inputVolume.get().getSamplingRate()
            if (self.noiseonlyinhalves is False):
                params += ' --noiseonlyinhalves'
        else:
            params += ' --sampling_rate %f' % self.inputVolumes.get().getSamplingRate()
        params += ' --number_frequencies %f' % Nfreqs
        params += ' --minRes %f' % self.minRes.get()
        params += ' --maxRes %f' % self.maxRes.get()
        params += ' --volumeRadius %f' % xdim
        params += ' --chimera_volume %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE_CHIMERA)
        params += ' --sym %s' % self.symmetry.get()
        params += ' --significance %f' % self.significance.get()
        params += ' --md_outputdata %s' % self._getExtraPath('mask_data.xmd')  
        if self.filterInput.get():
            params += ' --filtered_volume %s' % self._getExtraPath(FN_FILTERED_MAP)
        else:
            params += ' --filtered_volume %s' % ''

        self.runJob('xmipp_resolution_monogenic_signal', params)


    def createHistrogram(self):

        params = ' -i %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        params += ' --mask binary_file %s' % self._getExtraPath(OUTPUT_MASK_FILE)
        params += ' --steps %f' % 30
        params += ' --range %f %f' % (self.min_res_init, self.max_res_init)#(self.minRes.get(), self.maxRes.get())
        params += ' -o %s' % self._getExtraPath('hist.xmd')

        self.runJob('xmipp_image_histogram', params)
        
        
    def readMetaDataOutput(self):
        mData = md.MetaData(self._getExtraPath(METADATA_MASK_FILE))
        NvoxelsOriginalMask = float(mData.getValue(md.MDL_COUNT, mData.firstObject()))
        NvoxelsOutputMask = float(mData.getValue(md.MDL_COUNT2, mData.firstObject()))
        nvox = int(round(((NvoxelsOriginalMask-NvoxelsOutputMask)/NvoxelsOriginalMask)*100))
        return nvox

    def getMinMax(self, imageFile):
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        min_res = round(np.amin(imgData) * 100) / 100
        max_res = round(np.amax(imgData) * 100) / 100
        return min_res, max_res

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
            
        #Setting the min max for the summary
        imageFile = self._getExtraPath(OUTPUT_RESOLUTION_FILE_CHIMERA)
        min_, max_ = self.getMinMax(imageFile)
        self.min_res_init.set(round(min_*100)/100)
        self.max_res_init.set(round(max_*100)/100)
        self._store(self.min_res_init)
        self._store(self.max_res_init)

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

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'outputVolume'):
            messages.append(
                'Information about the method/article in ' + MONORES_METHOD_URL)
        return messages
    
    def _summary(self):
        summary = []
        summary.append("Highest resolution %.2f Å,   "
                       "Lowest resolution %.2f Å. \n" % (self.min_res_init,
                                                         self.max_res_init))
        Nvox = self.readMetaDataOutput()

        if (Nvox>10):
            summary.append("The resolution of %i %% of the mask voxels could not be computed. Maybe the mask was"
            "not correctly created, it is too wide or the resolution range does not cover the resolution at those voxels. "
            "If it is not the problem, decrease the significance in the advaced parameters can be an alternative" % Nvox)

        return summary

    def _citations(self):
        return ['Vilas2017']

