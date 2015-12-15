'''
Created on 5th Feb, 2014

@author: Roberto Marabini
         J.M. de la Rosa Trevin
'''
from pyworkflow.object import Float
from pyworkflow.tests import BaseTest, setupTestOutput, DataSet
from pyworkflow.em.data import SetOfParticles, CTFModel, Acquisition, Coordinate, Particle

import pyworkflow.em.metadata as md
from pyworkflow.em.packages.relion.convert import convertBinaryFiles
import pyworkflow.em.packages.relion as relion


class TestConversions(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')  
        cls.getFile = cls.dataset.getFile
        
        cls.ds = DataSet.getDataSet('relion_tutorial')
        
        
    def aaatest_particlesToStar(self):
        """ Write a SetOfParticles to Relion star input file. """
        imgSet = SetOfParticles(filename=self.getOutputPath("particles.sqlite"))
        n = 10
        fn = self.particles
        ctfs = [CTFModel(defocusU=10000, defocusV=15000, defocusAngle=15),
                CTFModel(defocusU=20000, defocusV=25000, defocusAngle=25)
               ]
        acquisition = Acquisition(magnification=60000, voltage=300,
                                  sphericalAberration=2., amplitudeContrast=0.07)
        imgSet.setAcquisition(acquisition)
        coord = Coordinate()
        coord.setMicId(1)

        for i in range(n):
            p = Particle()
            p.setLocation(i+1, fn)
            ctf = ctfs[i%2]
            p.setCTF(ctf)
            p.setAcquisition(acquisition)
            p._xmipp_zScore = Float(i)
            coord.setX(i*10)
            coord.setY(i*10)
            p.setCoordinate(coord)
            imgSet.append(p)
            
        fnStar = self.getOutputPath('particles.star')
        fnStk = self.getOutputPath('particles.stk')
        
        print ">>> Writing to file: %s" % fnStar
        relion.writeSetOfParticles(imgSet, fnStar, fnStk)
        
        mdAll = md.MetaData(fnStar)
        self.assertTrue(mdAll.containsLabel(md.RLN_IMAGE_COORD_X))
        self.assertTrue(mdAll.containsLabel(md.RLN_IMAGE_COORD_Y))
        self.assertFalse(mdAll.containsLabel(md.RLN_SELECT_PARTICLES_ZSCORE))
        
    def test_particlesFromStar(self):
        """ Read a set of particles from an .star file.  """
        fnStar = self.getFile('relion_it020_data')
        
        print ">>> Reading star file: ", fnStar
        mdAll = md.MetaData(fnStar)
        goldLabels = ['rlnVoltage', 'rlnDefocusU', 'rlnDefocusV', 
                      'rlnDefocusAngle', 'rlnSphericalAberration', 
                      'rlnAmplitudeContrast', 'rlnImageName', 'rlnImageId', 
                      'rlnCoordinateX', 'rlnCoordinateY', 'rlnMagnificationCorrection',
                      'rlnNormCorrection', 'rlnMicrographName', 'rlnGroupNumber', 
                      'rlnOriginX', 'rlnOriginY', 'rlnAngleRot', 'rlnAngleTilt', 
                      'rlnAnglePsi', 'rlnClassNumber', 'rlnLogLikeliContribution', 
                      'rlnNrOfSignificantSamples', 'rlnMaxValueProbDistribution']
        self.assertEqual(goldLabels, [md.label2Str(l) for l in mdAll.getActiveLabels()])
        self.assertEqual(4700, mdAll.size())
        

class TestConvertBinaryFiles(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.ds = DataSet.getDataSet('xmipp_tutorial')  
        cls.dsEmx = DataSet.getDataSet('emx')
        
    def test_hdfToStk(self):
        """ In this case the hdf stack files should be converted
        to .stk spider files for Relion.
        """
        stackFiles = ['BPV_1386_ptcls.hdf',
                      'BPV_1387_ptcls.hdf',
                      'BPV_1388_ptcls.hdf']
        
        partSet = SetOfParticles(filename=':memory:')
        
        for fn in stackFiles:
            particle = Particle()
            particle.setLocation(1, self.ds.getFile('particles/%s' % fn))
            partSet.append(particle)
            
        outputDir = self.getOutputPath()
        
        filesDict = convertBinaryFiles(partSet, outputDir)
        
        partSet.close()
        
        print filesDict
        
    def test_mrcsLink(self):
        """ In this case just a link with .mrcs extension 
        should be created
        """
        stackFile = self.dsEmx.getFile('particles/particles.mrc')
        partSet = SetOfParticles(filename=':memory:')
        
        for i in range(1, 10):
            particle = Particle()
            particle.setLocation(i, stackFile)
            partSet.append(particle)
            
        outputDir = self.getOutputPath()
        
        filesDict = convertBinaryFiles(partSet, outputDir)
        
        print filesDict        
    
            
        
        
        
