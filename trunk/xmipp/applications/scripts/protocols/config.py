sections = [
('Preprocessing', 
   [['Preprocess Micrograph', 'Preprocess Micrograph'], 
    ['Particles picking', 'Particles picking'], 
    ['Preprocess Particles', 'Preprocess Particles']]),
('2D', 
   [['Align+Classify', 'ML2D', 'CL2D'], 
    ['Align', 'ML2D', 'CL2D'], 
    ['Classify', 'KerDenSOM', 'Rotational Spectra']]),
('3D', 
   [['Initial Model', 'Common Lines', 'Random Conical Tilt'], 
    ['Model Refinement', 'Projection Matching']]),
('Other', [['Browse']])]

class ProtocolNames:
    preprocess_micrographs = 'preprocess_micrographs'    
    particle_pick = 'particle_pick' 
    preprocess_particles = 'preprocess_particles' 
    ml2d = 'ml2d'
    cl2d = 'cl2d'
    kerdensom = 'kerdensom'
    rotspectra = 'rotspectra'
    commonlines = 'commonlines'
    rct = 'rct'
    projmatch = 'projmatch'

launchDict = {
              'Preprocess Micrograph': ProtocolNames.preprocess_micrographs,
              'Particles picking':     ProtocolNames.particle_pick, 
              'Preprocess Particles':  ProtocolNames.preprocess_particles, 
              'ML2D':                  ProtocolNames.ml2d,
              'CL2D':                  ProtocolNames.cl2d,
              'KerDenSOM':             ProtocolNames.kerdensom,
              'Rotational Spectra':    ProtocolNames.rotspectra,
              'Common Lines':          ProtocolNames.commonlines,
              'Random Conical Tilt':   ProtocolNames.rct,
              'Projection Matching':   ProtocolNames.projmatch
              }
projectDefaults = {
                   'Cfg': '.project.cfg',
                   'Db': '.project.sqlite',
                   'LogsDir': 'Logs',
                   'RunsDir': 'Runs',
                   'RunsPrefix': 'run',
                   'TableGroups': 'groups',
                   'TableParams': 'params',
                   'TableProtocols': 'protocols',
                   'TableProtocolsGroups': 'protocols_groups',
                   'TableRuns': 'runs',
                   'TableSteps': 'steps',
                   'TableStepsRestart': 'steps_restart'
                   } 
    