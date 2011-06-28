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

launchDict = {
              'Preprocess Micrograph': 'preprocess_micrographs',
              'Particles picking':     'particle_pick', 
              'Preprocess Particles':  'preprocess_particles', 
              'ML2D':                  'ml2d',
              'CL2D':                  'cl2d',
              'KerDenSOM':             'kerdensom',
              'Rotational Spectra':    'rotspectra',
              'Common Lines':          'commonlines',
              'Random Conical Tilt':   'rct',
              'Projection Matching':   'projmatch'
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
    