# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Grigory Sharov (sharov@igbmc.fr)
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

import pyworkflow.em.packages.imagic as imagic
from pyworkflow.em.protocol import ProtImportParticles

from pyworkflow.tests import setupTestProject, DataSet
from test_workflow import TestWorkflow

  
   
class TestImagicWorkflow(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.particlesFn = cls.dataset.getFile('particles/xmipp_particles.xmd')
        
    def validateFilesExist(self, files):
        exist = []
        for f in files:
            if os.path.exists(f):
                exist.append(f)
        self.assertEqual(files, exist, "Missing output files")
    
    def test_msaWorkflow(self):
        """ Run an Import particles protocol. """
        protImport = self.newProtocol(ProtImportParticles,
                                      importFrom=2,
                                      mdFile=self.particlesFn,
                                      samplingRate=3.5)
        self.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. '
                            'outputParticles is None.' % self.particlesFn)

        protMsa = self.newProtocol(imagic.ImagicProtMSA,
                                      objLabel='imagic - msa',
                                      numberOfFactors=10,
                                      numberOfIterations=5,
                                      numberOfMpi=1)
        protMsa.inputParticles.set(protImport.outputParticles)
        self.launchProtocol(protMsa)


        protMsaClassify = self.newProtocol(imagic.ImagicProtMSAClassify,
                                          objLabel='imagic - msa classify',
                                          numberOfFactors=5,
                                          numberOfClasses=4)
        protMsaClassify.inputMSA.set(protMsa)
        self.launchProtocol(protMsaClassify)
        self.assertIsNotNone(protMsaClassify.outputClasses,
                             "There was a problem with the MSA-classify protocol's "
                             "outputClasses")