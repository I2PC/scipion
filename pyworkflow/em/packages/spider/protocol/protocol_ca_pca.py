# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from os.path import join

from pyworkflow.protocol.params import IntParam, PointerParam, EnumParam, FloatParam
from pyworkflow.em.convert import ImageHandler

from ..constants import CA
from ..spider import PcaFile
from protocol_base import SpiderProtocol



class SpiderProtCAPCA(SpiderProtocol):
    """Correspondence Analysis (CA) or Principal Component Analysis (PCA).
    
    CA is the preferred method of finding inter-image variations. 
    PCA computes the distance between data vectors with Euclidean 
    distances, while CA uses Chi-squared distance.     
    CA is superior because it ignores differences in exposure 
    between images, eliminating the need to rescale between images. 
    In contrast, PCA seems to be more robust: 
    less likely to be trapped in an infinite loop of numerical inaccuracy.
    
    For more info see:
    [[http://spider.wadsworth.org/spider_doc/spider/docs/techs/classification/tutorial.html#CAPCA][SPIDER MDA documentation]] 
    """
    _label = 'capca'
    
    def __init__(self, **kwargs):
        SpiderProtocol.__init__(self, **kwargs)
        self._caDir = 'CA'
        self._caPrefix = 'cas' 
        
        caFilePrefix = join(self._caDir, self._caPrefix + '_')
        
        self._params = {'ext': 'stk',
                        'particles': 'input_particles',
                        'particlesSel': 'input_particles_sel',
                        'outputParticles': 'output_particles',
                        'mask': 'mask',
                        # TO DO: read tags in case filenames change in SPIDER procedure
                        'imcFile': caFilePrefix + 'IMC',
                        'seqFile': caFilePrefix + 'SEQ',
                        'eigFile': caFilePrefix + 'EIG',
                        'eigenimages': join(self._caDir, 'stkeigenimg'),
                        'reconstituted': join(self._caDir, 'stkreconstituted')
                        }
    
    #--------------------------- DEFINE param functions --------------------------------------------  
     
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input particles to perform CA, PCA, or IPCA.')        
        form.addParam('analysisType', EnumParam, default=CA, choices=['CA', 'PCA', 'IPCA'],
                      label='Analysis type',
                      help='Select which type of analysis you want to perform: \n'
                           'correspondence analysis (CA), principal component analysis (PCA), '
                           'or iterative principal component analysis (IPCA)')
        form.addParam('addConstant', FloatParam, default=0,
                      condition="analysisType==%d" % CA, 
                      label='Additive constant',
                      help='Correspondence analysis requires the data to be positive. '
                           'In the case of negative values, a constant needs to be added. '
                           'An additive constant of *0* means automatic.')       
        form.addParam('numberOfFactors', IntParam, default=25,
                      label='Number of factors',
                      help='A 64x64 image can be expressed as a vector of 4096 dimensions. '
                           'In this step, we will reduce this number of dimensions to the number of factors specified here. '
                           'These factors will represent the largest systematic variations in the data.')
        form.addParam('maskType', EnumParam, 
                      choices=['circular', 'object'], default=0, 
                      display=EnumParam.DISPLAY_HLIST,
                      label='Mask type', 
                      help='Select which type of mask do you want to apply. '
                           'Only the pixels beneath this mask will be analyzed. '
                           'In the simplest case, a circular mask can be used. '
                           'Alternatively, a custom mask can be used '
                           'which follows the contour of the particle (but not too tightly).')
        form.addParam('radius', IntParam, default=-1,
                      label='Mask radius (px)', condition='maskType==0',
                      help='If -1, the entire image (in pixels) will be considered.')
        form.addParam('maskImage', PointerParam, label="Mask image", condition='maskType==1',
                      pointerClass='Mask', 
                      help="Select a mask file")       
        
    #--------------------------- INSERT steps functions --------------------------------------------  
    
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles',
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        if self.maskType > 0:
            self._insertFunctionStep('convertMaskStep', self.maskImage.get().getObjId())
        else:
            self.maskImage.set(None)
            
        self._insertFunctionStep('capcaStep', self.analysisType.get(), 
                                 self.numberOfFactors.get(), self.maskType.get())
        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions --------------------------------------------    
       
    def convertMaskStep(self, maskType):
        """ Convert the input mask if needed and
        copy some spider needed scripts. 
        """
        # Copy mask if selected
        if maskType > 0: # mask from file
            maskFn = self._getFileName('mask')
            ImageHandler().convert(self.maskImage.get().getLocation(), 
                       (1, maskFn))
        
    def capcaStep(self, analysisType, numberOfFactors, maskType):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        dim = self.inputParticles.get().getDimensions()[0]
        
        self._params.update({'[idim]': dim,
                             '[radius]': self.radius.get(),
                             '[cas-option]': analysisType + 1, # Index starts at 0
                             '[add-constant]': self.addConstant.get(),
                             '[num-factors]': numberOfFactors,
                             '[selection_doc]': self._params['particlesSel'],
                             '[particles]': self._params['particles'] + '@******',
                             '[custom_mask]': self._params['mask'] + '@1',
                             '[ca_dir]': self._caDir,
                             '[eigen_img]': self._params['eigenimages'], 
                             '[reconstituted_img]': self._params['reconstituted']
                             })
                   
        self.runTemplate('mda/ca-pca.msa', self.getExt(), self._params)
        
    def createOutputStep(self):
        # Generate outputs
        imc = PcaFile()
        imc.filename.set(self._getFileName('imcFile'))
        
        seq = PcaFile()
        seq.filename.set(self._getFileName('seqFile'))
        
        self._defineOutputs(imcFile=imc, seqFile=seq)        
        self._defineSourceRelation(self.inputParticles, imc)
        self._defineSourceRelation(self.inputParticles, seq)
        
            
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _summary(self):
        summary = []
        
        if self.analysisType == 0:
            summary.append(    'Analysis type: *Correspondence analysis*')
            if self.addConstant != 0:
                summary.append('    Additive constant: *%s*' % self.addConstant)
            else:
                summary.append('    Additive constant: *Auto*')
        if self.analysisType == 1:
            summary.append(    'Analysis type: *Principal component analysis*')
        if self.analysisType == 2:
            summary.append(    'Analysis type: *Iterative principal component analysis*')

        summary.append('Number of factors: *%s*' % self.numberOfFactors)
        
        if self.maskType == 0:  # circular mask
            if self.radius == -1:
                summary.append('Mask: *Circular, of radius 1/2 image dimension*')
            else:
                summary.append('Mask: *Circular, of radius: %s*' % self.radius)
        else:  # custom mask
            summary.append('Mask: *Custom file*')

        return summary
    
    def _methods(self):
        #print
        #print "protocol_ca_pca._methods.analysisType: %s" % self.analysisType 
        #print 
        
        msg  = "\nInput particles %s were subjected to " % self.getObjectTag('inputParticles')
        
        if self.analysisType == 0:
            msg += "correspondence analysis, "
        if self.analysisType == 1:
            msg += "principal component analysis, "
        if self.analysisType == 2:
            msg += "iterative principal component analysis, "
        
        msg += "computing %s factors, and using a " % self.numberOfFactors

        if self.maskType == 0:  # circular mask
            if self.radius == -1:
                msg += "circular mask of radius half the image dimension."
            else:
                msg += "circular mask of radius %s pixels." % self.radius
        else:  # custom mask
            msg += "custom mask %s." % self.getObjectTag('maskImage')
        
        return [msg]

        