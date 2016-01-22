#!/usr/bin/env xmipp_python
from xmipp import *
from protlib_base import *
import sys
import os

#init project
project = XmippProject()
#load project: read config file and open connection database
project.load()
# Clear project
project.clean()

pe = ProtocolExecutor('custom', project, script='custom_create_phantom_micrograph.py')
# Set values
valuesDict = {
  'Comment': 'Automatic test for micrograph generation',
  'AlphaUListStr': "10 20 30", 'AlphaTListStr': "40 50 60", 'TiltAngle': 65, 
  'NumberOfElements':73, 'NumberOfMicrographs': 3, 
  'Coordinates': 'random', 'CoordinatesFile': 'coords.txt',
  'DoFlip': True, 'FlipRate': 0.4, 'NoiseSigma': 3., 'DoRotate':True,
  'AddNoiseU': True, 'AddNoiseT': True
}
pe.setValues(valuesDict)
pe.runProtocol()

pe = ProtocolExecutor('import_micrographs', project)
# Set values
valuesDict = {
  'DirMicrographs': "./Custom/create_phantom_micrograph_001", 'ExtMicrographs': '*.mrc',
  'TiltPairs': True, 'PairDescr': 'tilted_pairs.xmd', 
  'SamplingRate': '1', 'NumberOfMpi': 1,   
}
pe.setValues(valuesDict)
pe.runProtocol()


pe = ProtocolExecutor('particle_pick', project, script='template/particle_pick_run_001.py')
pe.setValue('LaunchGUI', False)
pe.runProtocol()

os.system('cp Custom/create_phantom_micrograph_001/micrograph???[UT].pos ParticlePicking/Manual/run_001/extra/')
os.system('cp template/families.xmd ParticlePicking/Manual/run_001/extra/')
os.system('cp tilted_pairs.xmd ParticlePicking/Manual/run_001/')
os.system('sed -i "s/LaunchGUI = False/LaunchGUI = True/g"  Runs/particle_pick_run_001.py')

pe = ProtocolExecutor('extract_particles', project)
pe.setValues( {'Family': 'DefaultFamily', 'ParticleSize': 50, 
'DoFlip': False, 'DoRemoveDust': False,
'PickingRun': 'particle_pick_run_001',
'NumberOfMpi': 2, })
pe.runProtocol()

#pe = ProtocolExecutor('cl2d', project, script='template/cl2d_run_001.py')
#pe.setValue('
#pe.runProtocol()

pe = ProtocolExecutor('rct', project)
pe.setValues({'ExtractionRun': "extract_particles_run_001", 
'ClassifMd': 'classes.xmd', 'CenterMaxShift': 30,
'SelectedClasses': "1", 'NumberOfMpi': 2})
pe.runProtocol()
