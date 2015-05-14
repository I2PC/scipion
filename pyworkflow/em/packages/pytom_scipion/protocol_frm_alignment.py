# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin  (jmdelarosa@cnb.csic.es)
# *              Yuxiang Chen 
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as em  

from convert import printLicense, writeVolumesSqlite, writeSetOfVolumes



class ProtFrmAlign(em.ProtAlignVolume):
    """ Subtomogram averaging using pytom autofocus """
    _label = 'frm 3d align'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineParams(self, form):
        
        form.addSection('Input')
        form.addParam('inputVolumes', params.PointerParam, 
                      pointerClass='SetOfVolumes',
                      label="Input volume particles", important=True, 
                      help='Subtomograms to average')
        form.addParam('provideReference', params.BooleanParam, default=True,
                      label='Provide reference?',
                      help='Use _Yes_ if you want to provide you own\n'
                           'set of volumes as initial references.')
        form.addParam('inputReference', params.PointerParam,
                      pointerClass='Volume', allowsNull=True,
                      condition='provideReference',
                      label='Input reference volume',
                      help='Volume used as reference for alignment')
        form.addParam('referenceGenerate', params.LabelParam,
                       condition='not provideReference',
                       label='Initial reference will be generated with a random angular assigment'),
        form.addParam('numberOfIterations', params.IntParam, default=10,
                      label='Number of iterations',
                      help='')
        form.addParam('binning', params.IntParam, default=1, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Binning factor (int)')
        
        group = form.addGroup('Alignment')
        group.addParam('maxFreq', params.FloatParam, default=0.15,
                      label='Max. frequency (<0.5)',
                      help='the starting frequency (in digital frequency normalized to 0.5) at which the \n'
                           'alignment procedure would start')
        group.addParam('alignmentMask', params.PointerParam,
                       pointerClass='VolumeMask', allowsNull=True,
                       label='Alignment mask',
                       help='Mask used during alignment')
        group.addParam('offset', params.IntParam, 
                      label='Alignment offset',
                      help='The maximal spatial range (radius, in pixel) that   \n'
                           'the subtomogram would be shifted. This is necessary \n'
                           'to prevent shifting the volume out-of-frame and     \n'
                           'reduce the search space.')
        
        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        """ Mainly prepare the command line for calling reconstruct_significant program"""
        self.volXml = self._getExtraPath('volume_particles.xml')
        self.refVol = self._getExtraPath('input_reference.mrc')
        self.maskVol = self._getExtraPath('input_mask.mrc')
        
        self._insertFunctionStep('convertInputStep')
        self._insertFrmAlignStep()
        self._insertFunctionStep('createOutputStep') 
        
    def _insertFrmAlignStep(self):
        """ Prepare the dictionary to be used in job.xml template. """
        
        Ts=self.inputVolumes.get().getSamplingRate()
        xdim=self.inputVolumes.get().getDim()[0]
        args = {'particleList': os.path.basename(self.volXml),
                'reference': os.path.basename(self.refVol),
                'mask': os.path.basename(self.maskVol),
                'frequency': int(self.maxFreq.get()*xdim),
                'offset': self.offset.get(),
                'pixelSize': Ts,
                'diameter': xdim*Ts
                }
                
        self._insertFunctionStep('runFrmAlignStep', args)
        

    #--------------------------- STEPS functions --------------------------------------------        
    def convertInputStep(self):
        printLicense()
        self.info('Writing pytom xml file: ' + self.volXml)
        volsDir = self._getExtraPath('inputVolumes')
        pwutils.makePath(volsDir)
        writeSetOfVolumes(self.inputVolumes.get(), self.volXml, volsDir)
        ih = em.ImageHandler()
        
        self.info('Converting input reference to: ' + self.refVol)
        ih.convert(self.inputReference.get(), self.refVol)
        
        self.info('Converting input mask to: ' + self.maskVol)
        ih.convert(self.alignmentMask.get(), self.maskVol)
            
    def runFrmAlignStep(self, args):
        """ Prepare the job.xml file and run FRMAlignment.py script. """
        jobTemplate = """
<FRMJob Destination='.' BandwidthRange='[4, 64]' Frequency='%(frequency)s' MaxIterations='10' PeakOffset='10' AdaptiveResolution='0.1' FSC='0.5'>

<Reference PreWedge="" File="%(reference)s" Weighting="">
</Reference>
<Mask Filename="%(mask)s" Binning="1"/>
<SampleInformation PixelSize="%(pixelSize)s" ParticleDiameter="%(diameter)s"/>

<ParticleListLocation Path="%(particleList)s"/>

</FRMJob> 
        """
#<Mask Filename="%(mask)s" Binning="1" isSphere="True"/>
        jobName = 'job.xml'
        jobXml = self._getExtraPath(jobName)
        self.info('Writing job.xml file to: ' + jobXml)
        f = open(jobXml, 'w')
        f.write(jobTemplate % args)
        f.close()
        script = self._getScript("frm", "FRMAlignment.py")
        # Run the alignment with the job script (-j option) 
        # and in verbose mode (-v option)
        self.runJob('python', '%s -j %s -v ' % (script, jobName), 
                    cwd=self._getExtraPath())
        #print('python', '%s -j %s' % (script, jobName))
        # scipion run python /home/coss/usb_linux/scipion/scipion/software/em/pytom/classification/auto_focus_classify.py -p volumes.xml -k 1 -f 2 -m ../phantom.mrc -s 5 -i 10 -n 0.2 -g -2 -t 0.4

    def createOutputStep(self):
        pass

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        """ Check that some preconditions are met before launching 
        the auto-focus classification run. 
        """
        errors = []
        
        if self.numberOfMpi < 0:
            errors.append('Number of MPI should be greater than 2.')
            
        if self.provideReference:
            if self.inputReference.get() is None:
                errors.append('Please select the input reference.')
        else:
            errors.append('Not implemented!!!')
            
            
        inputVols = self.inputVolumes.get()
        if inputVols is not None:
            xdim = inputVols.getDim()[0]
            half = xdim/2
            
            if self.maxFreq >= half:
                errors.append('Frequency should be less dim/2 pixels (%d/2=%d)' % (xdim, half))
        
        return errors
        
    def _summary(self):
        summary = []
        #summary.append("Input classes: %s" % self.inputReferences.get().getNameId())
        #summary.append("Starting from: %d random volumes" % self.numberOfReferences)
        return summary
    
    def _citations(self):
        return ['Chen2014']
    
    def _methods(self):
        return []
    
    
    #--------------------------- UTILS functions --------------------------------------------   
    
    def _getScript(self, *paths):
        return os.path.join(os.environ['PYTOM_HOME'], *paths)
    
    def getAverageMap(self, it):
        return self._getExtraPath('average_iter%d.em' % it)
    
    def getVolumesSqlite(self, it):
        volsXml = self._getExtraPath('aligned_pl_iter%d.xml' % it)
        volsSqlite = volsXml + '.sqlite'
        
        if not os.path.exists(volsSqlite):
            writeVolumesSqlite(volsXml, volsSqlite)
            
        return volsSqlite
    
