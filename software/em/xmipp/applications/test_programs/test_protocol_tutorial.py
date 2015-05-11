#!/usr/bin/env xmipp_python
import os, shutil, sys
from os.path import exists, join, abspath
from datetime import datetime, timedelta
from protlib_xmipp import greenStr, warnStr, redStr, XmippScript
from protlib_base import ProtocolExecutor, XmippProject


class ScriptProtocolTester(XmippScript):    
    def __init__(self):
        XmippScript.__init__(self,True)
        
    def defineParams(self):
        self.addUsageLine('Run program tests')
        ## params
        self.addParamsLine(' [-d <directory=tmpProject>]   : OutputDirectory')
        self.addParamsLine('   alias --output;')

    def run(self):
        import urllib
        fn = "xmipp3_intro.tgz"
        
        path = self.getParam("-d")
        
        if not exists(path):
            os.makedirs(path)
        # Change to output directory
        os.chdir(path)
        
        url = "http://xmipp.cnb.csic.es/Downloads/XmippDoc/3.1/%s" % fn
        if not exists(fn):
            print "Downloading from: %s" % url
            urllib.urlretrieve(url, fn)
        
        fnProject = "BPV_Project"
        if not exists(fnProject):
            print "Extracting project folder"
            os.system('tar -xvzf %s' % fn)
        
        os.chdir(fnProject)
        
        #Create project
        project = XmippProject()
        #project.load()
        project.clean()
        
        pe = ProtocolExecutor('import_micrographs', project)
        pe.setValues({'DirMicrographs': "InputData", 'ExtMicrographs': '*.mrc',
                      'Voltage': 300, 'SphericalAberration': 1.2,
		      'SamplingRate': '1.237', 'NumberOfMpi': 1, 'SubmitToQueue': 0})
        pe.runProtocol()
        
        pe = ProtocolExecutor('screen_micrographs', project)
        pe.setValues({'ImportRun': 'import_micrographs_run_001',
		'DownsampleFactor': 2, 'NumberOfMpi': 4, 'SubmitToQueue': 0})
        pe.runProtocol()   
        
        pe = ProtocolExecutor('downsample_micrographs', project)
        pe.setValues({'ImportRun': 'screen_micrographs_run_001',
		'DownsampleFactor': 5, 'NumberOfMpi': 4, 'SubmitToQueue': 0})
        pe.runProtocol()       
        
        
        pe = ProtocolExecutor('particle_pick', project)
        pe.setValues({'ImportRun': 'downsample_micrographs_run_001',
                      'LaunchGUI': False})
        pe.runProtocol()
        # Copy coordinates and family metadata   
        os.system('cp ../../coordinates/* ParticlePicking/Supervised/run_001/extra/')    
        # Reset LaunchGUI to True 
        os.system('sed -i "s/LaunchGUI = False/LaunchGUI = True/g"  Runs/particle_pick_run_001.py') 
        
        pe = ProtocolExecutor('extract_particles', project)
        pe.setValues({'PickingRun': 'particle_pick_run_001',
                      'ParticleSize': 110,
		      'DoInvert': True, 'NumberOfMpi': 4, 'SubmitToQueue': 0})
        pe.runProtocol() 
        
        pe = ProtocolExecutor('projmatch', project)
        pe.setValues({'SelFileName': "Images/Extracted/run_001/images.xmd",
                      'ReferenceFileNames': "../../InitialVolume/BPV_scale_filtered_windowed.vol",
                      'MaskRadius': 55, 'InnerRadius': "28", 'OuterRadius': "55",
                      'SymmetryGroup': "i1", 
		      'MpiJobSize': "5", 'NumberOfMpi': 4, 'SubmitToQueue': 0})
        pe.runProtocol() 

if __name__ == '__main__':
    ScriptProtocolTester().tryRun()

