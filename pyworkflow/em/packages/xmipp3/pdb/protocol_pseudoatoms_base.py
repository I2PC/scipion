# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
# *           Slavica Jonic                (jonic@impmc.upmc.fr)
# * Ported to Scipion:
# *           J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Jan 2014
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

from os.path import basename

from pyworkflow.object import String
from pyworkflow.utils.path import cleanPattern, createLink, moveFile
from pyworkflow.protocol.params import EnumParam, PointerParam, FloatParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.protocol import Prot3D
from pyworkflow.em.packages.xmipp3.convert import getImageLocation

#from xmipp3 import XmippProtocol
NMA_MASK_NONE = 0
NMA_MASK_THRE = 1
NMA_MASK_FILE = 2



class XmippProtConvertToPseudoAtomsBase(Prot3D):
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addParam('maskMode', EnumParam, choices=['none', 'threshold', 'file'], 
                      default=NMA_MASK_NONE, 
                      label='Mask mode', display=EnumParam.DISPLAY_COMBO,
                      help='')        
        form.addParam('maskThreshold', FloatParam, default=0.01, 
                      condition='maskMode==%d' % NMA_MASK_THRE,
                      label='Threshold value',
                      help='Gray values below this threshold are set to 0')
        form.addParam('volumeMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask volume', condition='maskMode==%d' % NMA_MASK_FILE,
                      )          
        form.addParam('pseudoAtomRadius', FloatParam, default=1, 
                      label='Pseudoatom radius (vox)',
                      help='Pseudoatoms are defined as Gaussians whose \n'
                           'standard deviation is this value in voxels') 
        form.addParam('pseudoAtomTarget', FloatParam, default=5,
                      expertLevel=LEVEL_ADVANCED, 
                      label='Volume approximation error(%)',
                      help='This value is a percentage (between 0.001 and 100) \n'
                           'specifying how fine you want to approximate the EM \n'
                           'volume by the pseudoatomic structure. Lower values \n'
                           'imply lower approximation error, and consequently, \n'
                           'more pseudoatoms.')        
             
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertMaskStep(self, fnVol, prefix=''):
        """ Check the mask selected and insert the necessary steps.
        Return the mask filename if needed.
        """
        fnMask = ''
        if self.maskMode == NMA_MASK_THRE:
            fnMask = self._getExtraPath('mask%s.vol' % prefix)
            maskParams = '-i %s -o %s --select below %f --substitute binarize' % (fnVol, fnMask, self.maskThreshold.get())
            self._insertRunJobStep('xmipp_transform_threshold', maskParams)
        elif self.maskMode == NMA_MASK_FILE:
            fnMask = getImageLocation(self.volumeMask.get())
        return fnMask
        
        
    #--------------------------- STEPS functions --------------------------------------------
    def convertToPseudoAtomsStep(self, inputFn, fnMask, sampling, prefix=''):
        pseudoatoms = 'pseudoatoms%s'%prefix
        outputFn = self._getPath(pseudoatoms)
        sigma = sampling * self.pseudoAtomRadius.get() 
        targetErr = self.pseudoAtomTarget.get()
        nthreads = self.numberOfThreads.get()
        params = "-i %(inputFn)s -o %(outputFn)s --sigma %(sigma)f --thr %(nthreads)d "
        params += "--targetError %(targetErr)f --sampling_rate %(sampling)f -v 2 --intensityColumn Bfactor"
        if fnMask:
            params += " --mask binary_file %(fnMask)s"
        self.runJob("xmipp_volume_to_pseudoatoms", params % locals())
        for suffix in ["_approximation.vol", "_distance.hist"]:
            moveFile(self._getPath(pseudoatoms+suffix), self._getExtraPath(pseudoatoms+suffix))
        cleanPattern(self._getPath(pseudoatoms+'_*'))
        
    def createChimeraScript(self, volume, pdb):
        """ Create a chimera script to visualize a pseudoatoms pdb
        obteined from a given EM 3d volume.
        A property will be set in the pdb object to 
        store the location of the script.
        """
        pseudoatoms = pdb.getFileName()
        scriptFile = pseudoatoms + '_chimera.cmd'
        pdb._chimeraScript = String(scriptFile)
        sampling = volume.getSamplingRate()
        radius = sampling * self.pseudoAtomRadius.get() 
        fnIn = volume.getFileName()
        localInputFn = self._getBasePath(fnIn)
        createLink(fnIn, localInputFn)
        fhCmd = open(scriptFile, 'w')
        fhCmd.write("open %s\n" % basename(pseudoatoms))
        fhCmd.write("rangecol bfactor,a 0 white 1 red\n")
        fhCmd.write("setattr a radius %f\n" % radius)
        fhCmd.write("represent sphere\n")
        fhCmd.write("open %s\n" % basename(localInputFn))
         
        threshold = 0.01
        if self.maskMode == NMA_MASK_THRE:
            self.maskThreshold.get()
        xdim = volume.getDim()[0]
        origin = xdim / 2
        fhCmd.write("volume #1 level %f transparency 0.5 voxelSize %f originIndex %d\n" % (threshold, sampling, origin))
        fhCmd.close()
     
