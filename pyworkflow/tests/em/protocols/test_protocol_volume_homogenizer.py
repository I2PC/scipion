# ***************************************************************************
# * Authors:     Mohsen Kazemi (mkazemi@cnb.csic.es)
# *              Javier Vargas (javier.vargasbalbuena@mcgill.ca)  
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
# ***************************************************************************/

import os
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportParticles, ProtImportVolumes, ProtSubSet 
from pyworkflow.em.packages.xmipp3.protocol_volume_homogenizer import XmippProtVolumeHomogenizer
from pyworkflow.em.packages.relion import ProtRelionRefine3D
from pyworkflow.em.data import SetOfParticles
    
class TestVolumeHomogenizer(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
           
    def importFromRelionRefine3D(self):
        """ Import aligned Particles
        """
        prot = self.newProtocol(ProtImportParticles,
                                objLabel='particles from relion (auto-refine 3d)',
                                importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                starFile=self.dsRelion.getFile('import/classify3d/extra/relion_it015_data.star'),
                                magnification=10000,
                                samplingRate=7.08,
                                haveDataBeenPhaseFlipped=True
                                )
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputParticles.getFileName(),
                             "There was a problem with the import")
    
        return prot
    
    def getParticles(self, protImport, classid):
        dbPartSet = protImport._getPath("particles_class-%d.sqlite" % classid)
        class3D = protImport.outputClasses[classid]
        if os.path.exists(dbPartSet):
            os.remove(dbPartSet)
        partSet = SetOfParticles(filename=dbPartSet)
        partSet.copyInfo(class3D)
    
        for part in class3D:
            partSet.append(part)
        partSet.write()
        partSet.close()
        
        protImportCls1 = self.newProtocol(ProtImportParticles,
                                objLabel='particles class-%d' % classid,
                                importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                sqliteFile=dbPartSet,
                                magnification=10000,
                                samplingRate=7.08,
                                haveDataBeenPhaseFlipped=True
                                )
        self.launchProtocol(protImportCls1)
        self.assertIsNotNone(protImportCls1.outputParticles.getFileName(),
                             "There was a problem with the import")
        return protImportCls1
        
    def subSetOfParticles(self,importPartClass1):
    
        protSubset = self.newProtocol(ProtSubSet,
                                       objLabel='subset of particles',
                                       chooseAtRandom=True,
                                       nElements=200)
    
        protSubset.inputFullSet.set(importPartClass1.outputParticles)
        self.launchProtocol(protSubset)
        self.assertIsNotNone(protSubset.outputParticles.getFileName(),
                              "There was a problem with the protSubset")
    
        return protSubset
    
 
    def getVolume(self,  classid):
        """ Import a Volumes
        """                        
        protImportVol = self.newProtocol(ProtImportVolumes,
                                objLabel='vol class-%d' % classid,
                                filesPath=self.dsRelion.getFile('import/classify3d/extra/relion_it015_class00%d.mrc' % classid),
                                samplingRate=7.08,
                                )
        self.launchProtocol(protImportVol)        
        return protImportVol
    
    def relionRefine(self,protSubSetClass, vol, classid):
            
            relionRefine = self.newProtocol(ProtRelionRefine3D,
                                            objLabel='relion_%s' % classid,
                                            doCTF=False, runMode=1,
                                            memoryPreThreads=1,
                                            maskDiameterA=424,
                                            angularSamplingDeg=1,
                                            localSearchAutoSamplingDeg = 2,
                                            symmetryGroup="c1",
                                            numberOfMpi=3, numberOfThreads=1)
            
            relionRefine.inputParticles.set(protSubSetClass.outputParticles)
            relionRefine.referenceVolume.set(vol.outputVolume)
                        
            self.launchProtocol(relionRefine)
            return relionRefine
     
    def homoNoGoldStandardNoAlign(self,relion1,relion2):
        '''test without Goldstandard and alignment'''
        protVolumeHomogenizer = self.newProtocol(XmippProtVolumeHomogenizer,
                                objLabel='volume homogenizer1')
        
        protVolumeHomogenizer.referenceVolume.set(relion1.outputVolume)
        protVolumeHomogenizer.inputVolume.set(relion2.outputVolume)
        protVolumeHomogenizer.inputParticles.set(relion2.outputParticles)
        protVolumeHomogenizer.doAlignment.set(False)  
        self.launchProtocol(protVolumeHomogenizer)
        
        self.assertIsNotNone(protVolumeHomogenizer.outputParticles.getFileName(),
                              "There was a problem with the homoNoGoldStandardNoAlign")
        
        return protVolumeHomogenizer
    
    
    def homoNoGoldStandardAlign(self,relion1,relion2):
        '''test without Goldstandard and alignment'''
        protVolumeHomogenizer = self.newProtocol(XmippProtVolumeHomogenizer,
                                objLabel='volume homogenizer2')
        
        protVolumeHomogenizer.referenceVolume.set(relion1.outputVolume)
        protVolumeHomogenizer.inputVolume.set(relion2.outputVolume)
        protVolumeHomogenizer.inputParticles.set(relion2.outputParticles)
        protVolumeHomogenizer.doAlignment.set(True)  
        self.launchProtocol(protVolumeHomogenizer)
        
        self.assertIsNotNone(protVolumeHomogenizer.outputParticles.getFileName(),
                              "There was a problem with the homoNoGoldStandardAlign")
        
        return protVolumeHomogenizer
    
    def homoGoldStandardNoAlign(self,relion1,relion2):
        '''test without Goldstandard and alignment'''
        protVolumeHomogenizer = self.newProtocol(XmippProtVolumeHomogenizer,
                                objLabel='volume homogenizer3',
                                doGoldStandard = True)
                
        protVolumeHomogenizer.referenceVolume1.set(relion1.outputVolume)
        protVolumeHomogenizer.inputVolume1.set(relion2.outputVolume)
        protVolumeHomogenizer.referenceVolume2.set(relion1.outputVolume)
        protVolumeHomogenizer.inputVolume2.set(relion2.outputVolume)
        protVolumeHomogenizer.inputParticles.set(relion2.outputParticles)
        protVolumeHomogenizer.doAlignment.set(False)
          
        self.launchProtocol(protVolumeHomogenizer)        
        self.assertIsNotNone(protVolumeHomogenizer.outputParticles01.getFileName(),
                              "There was a problem with the homoGoldStandardNoAlign")
        
        self.assertIsNotNone(protVolumeHomogenizer.outputParticles02.getFileName(),
                              "There was a problem with the homoGoldStandardNoAlign")
        
        return protVolumeHomogenizer
    
    def homoGoldStandardAlign(self,relion1,relion2):
        '''test without Goldstandard with alignment'''
        protVolumeHomogenizer = self.newProtocol(XmippProtVolumeHomogenizer,
                                objLabel='volume homogenizer4',
                                doGoldStandard = True)
                
        protVolumeHomogenizer.referenceVolume1.set(relion1.outputVolume)
        protVolumeHomogenizer.inputVolume1.set(relion2.outputVolume)
        protVolumeHomogenizer.referenceVolume2.set(relion1.outputVolume)
        protVolumeHomogenizer.inputVolume2.set(relion2.outputVolume)
        protVolumeHomogenizer.inputParticles.set(relion2.outputParticles)
        protVolumeHomogenizer.doAlignment.set(True)
          
        self.launchProtocol(protVolumeHomogenizer)        
        self.assertIsNotNone(protVolumeHomogenizer.outputParticles01.getFileName(),
                              "There was a problem with the homoGoldStandardAlign")
        
        self.assertIsNotNone(protVolumeHomogenizer.outputParticles02.getFileName(),
                              "There was a problem with the homoGoldStandardAlign")
        
        return protVolumeHomogenizer
            
        
    def test_volumeHomogenizer(self):
        """ Test protocol volume homogenizer
        """
        protImpClasses = self.importFromRelionRefine3D()
        
        importPartClass1 = self.getParticles(protImpClasses, classid=1)
        protSubSetClass1 = self.subSetOfParticles(importPartClass1)
    
        importPartClass2 = self.getParticles(protImpClasses, classid=2)
        protSubSetClass2 = self.subSetOfParticles(importPartClass2)
        
        vol1 = self.getVolume(classid=1)
        vol2 = self.getVolume(classid=2)
        
        relion1 = self.relionRefine(protSubSetClass1,vol1,classid=1)
        relion2 = self.relionRefine(protSubSetClass2,vol2,classid=2)
                
        homo1 = self.homoNoGoldStandardNoAlign(relion1,relion2)
        #homo2 = self.homoNoGoldStandardAlign(relion1,relion2)
        homo3 = self.homoGoldStandardNoAlign(relion1,relion2)
        #homo4 = self.homoGoldStandardAlign(relion1,relion2)


        

        
        
        
        