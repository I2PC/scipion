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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.protocol.params import (PointerParam, FloatParam,  
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume 
from pyworkflow.em.protocol import ProtReconstruct3D
from pyworkflow.em.packages.relion.convert import convertBinaryFiles

class ProtRelionReconstruct(ProtReconstruct3D):
    """    
    Reconstruct a volume using Relion from a given set of particles.
    The alignment parameters will be converted to a Relion star file
    and used as direction projections to reconstruct.
    """
    _label = 'reconstruct'
    ##doContinue = False
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      label="Input particles",  
                      help='Select the input images from the project.')     
#         form.addParam('doNormalize', BooleanParam, default=False,
#                       label='Normalize',
#                       help='If set to True, particles will be normalized in the way RELION prefers it.')
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Relion Symmetry][http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Relion') 
        form.addParam('maxRes', FloatParam, default=-1,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).\n'
                           'Param *--maxres* in Relion.') 
        form.addParam('pad', FloatParam, default=2,
                      label="Padding factor",  
                      help='Padding factor, *--pad* param in Relion.') 
        
        form.addParam('extraParams', StringParam, default='', expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters: ', 
                      help='Extra parameters to *relion_reconstruct* program:\n'
                      """
                      --subtract () : Subtract projections of this map from the images used for reconstruction                                                                                                                               
                       --NN (false) : Use nearest-neighbour instead of linear interpolation before gridding correction                                                                                                                       
                     --blob_r (1.9) : Radius of blob for gridding interpolation                                                                                                                                                              
                       --blob_m (0) : Order of blob for gridding interpolation                                                                                                                                                               
                      --blob_a (15) : Alpha-value of blob for gridding interpolation                                                                                                                                                         
                        --iter (10) : Number of gridding-correction iterations                                                                                                                                                               
                       --refdim (3) : Dimension of the reconstruction (2D or 3D)                                                                                                                                                             
               --angular_error (0.) : Apply random deviations with this standard deviation (in degrees) to each of the 3 Euler angles                                                                                                        
                 --shift_error (0.) : Apply random deviations with this standard deviation (in pixels) to each of the 2 translations
            --fom_weighting (false) : Weight particles according to their figure-of-merit (_rlnParticleFigureOfMerit)
                           --fsc () : FSC-curve for regularized reconstruction
                      """)
        form.addSection('CTF')#, condition='doCTF')
        form.addParam('doCTF', BooleanParam, default=False,
                       label='Apply CTF correction?',
                       help='Param *--ctf* in Relion.')
        form.addParam('ctfIntactFirstPeak', BooleanParam, default=False,
                      condition='doCTF',
                      label='Leave CTFs intact until first peak?',
                      help='Param *--ctf_intact_first_peak* in Relion.')
        form.addParam('onlyFlipPhases', BooleanParam, default=False,
                      condition='doCTF',
                      label='Do not correct CTF-amplitudes? (only flip phases)',
                      help='Param *--only_flip_phases* in Relion.')
        line = form.addLine('Beam tilt in direction: ',
                            condition='doCTF',
                            help='Beamtilt in the directions X and Y (in mrad)')
        line.addParam('beamTiltX', FloatParam, default='0.0', label='X ')
        line.addParam('beamTiltY', FloatParam, default='0.0', label='Y ')            
        
        form.addParallelSection(threads=2, mpi=0) 
        #TODO: Add an option to allow the user to decide if copy binary files or not        
            
    #--------------------------- INSERT steps functions --------------------------------------------  

    def _insertAllSteps(self): 
        ##self._initialize()
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertReconstructStep()
        self._insertFunctionStep('createOutputStep')

    def _getProgram(self, program='relion_refine'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    def _insertReconstructStep(self):
        imgSet = self.inputParticles.get()
        imgStar = self._getFileName('input_particles.star')
        
        params = ' --i %s' % imgStar
        params += ' --o %s' % self._getPath('output_volume.vol')
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --angpix %0.3f' % imgSet.getSamplingRate()
        params += ' --maxres %0.3f' % self.maxRes.get()
        params += ' --pad %0.3f' % self.pad.get()
        params += ' --j %d' % self.numberOfThreads.get()
        
        #TODO Test that the CTF part is working
        if self.doCTF:
            params += ' --ctf'
            if self.ctfIntactFirstPeak:
                params += ' --ctf_intact_first_peak'
            if self.onlyFlipPhases:
                params += ' --only_flip_phases' 
            params += ' --beamtilt_x %0.3f' % self.beamTiltX.get()
            params += ' --beamtilt_y %0.3f' % self.beamTiltY.get()
        
        self._insertFunctionStep('reconstructStep', params)

    #--------------------------- STEPS functions --------------------------------------------
    def reconstructStep(self, params):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        self.runJob(self._getProgram('relion_reconstruct'), params)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_particles.star': self._getTmpPath('input_particles.star'),
            'output_volume': self._getPath('output_volume.vol')
            }
        self._updateFilenamesDict(myDict)


    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet  = self.inputParticles.get()
        imgStar = self._getFileName('input_particles.star')

        from convert import writeSetOfParticles
        print "Before filesMapping"
        filesMapping = convertBinaryFiles(imgSet, self._getTmpDir())
        # Pass stack file as None to avoid write the images files
#        writeSetOfParticles(imgSet,imgStar,filesMapping,
#                            is2D=False, isInverseTransform=True,
#                            writeAlignment=True)
        writeSetOfParticles(imgSet,imgStar,filesMapping,
                            is2D=False, isInverseTransform=False,
                            writeAlignment=True)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        volume = Volume()
        volume.setFileName(self._getFileName('output_volume'))
        volume.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(imgSet, volume)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
