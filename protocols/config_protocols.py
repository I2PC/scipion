
#--------------------------------------------------------------------------------
# Protocols info (name, title, path)
#--------------------------------------------------------------------------------
protocols = {
        'import_micrographs': ('Import', 'Micrographs/Imported'),
        'downsample_micrographs': ('Downsample', 'Micrographs/Downsampled'),
        'screen_micrographs': ('Screen', 'Micrographs/Screen'),
        'particle_pick': ('Manual/Supervised', 'ParticlePicking/Supervised'),
        'particle_pick_auto': ('Automatic', 'ParticlePicking/Auto'),
        'extract_particles': ('Extract', 'Images/Extracted'),
        'import_particles': ('Import', 'Images/Imported'),
        'merge_particles': ('Merge', 'Images/Merged'),
        'preprocess_particles': ('Preprocess', 'Images/Preprocessed'),
        'screen_particles': ('Screen', 'Images/Screening'),
        'ml2d': ('ML2D', '2D/ML2D'),
        'cl2d': ('CL2D', '2D/CL2D'),
        'cl2d_align': ('Only align', '2D/Alignment'),
        'kerdensom': ('KerDenSOM',  '2D/KerDenSOM'),
        'rotspectra': ('Rotational Spectra', '2D/RotSpectra'),
        'screen_classes': ('Screen classes', '2D/Screening'),
        'rct': ('Random Conical Tilt', '3D/InitialVolume/RCT'),
        'initvolume_ransac': ('RANSAC', '3D/InitialVolume/RANSAC'),
        'preprocess_volume': ('Preprocess', '3D/InitialVolume/Preprocessed'),
        'create_volume_mask': ('Create Volume mask', '3D/Mask'),
        'projmatch': ('Projection Matching', '3D/ProjMatch'), 
        'ml3d': ('ML3D', '3D/ML3D'),
        'nma': ('Normal Mode Analysis', '3D/NMA'),
        'nma_alignment': ('Flexible alignment', '3D/NMA_alignment'),
        'relion3d': ('Relion3D', '3D/Relion3D'),
        'mltomo': ('MLTomo', '3D/MLTomo'),
        'subtraction': ('Partial Projection Subtraction', '3D/ProjSub'),
        'custom': ('Custom', 'Custom'),
        'xmipp': ('Xmipp Programs', 'XmippPrograms')            
        }

#--------------------------------------------------------------------------------
# Protocols sections and groups
#--------------------------------------------------------------------------------
sections = [
('Preprocessing', 
   [['Micrographs', 'import_micrographs','screen_micrographs','downsample_micrographs'], 
    ['Particle picking', 'particle_pick', 'particle_pick_auto'], 
    ['Particles', 'extract_particles', 'import_particles', 'merge_particles', ['Other', 'preprocess_particles', 'screen_particles']]]),
('2D', 
   [['Align+Classify', 'cl2d', 'ml2d', ['Other', 'cl2d_align', 'kerdensom', 'rotspectra', 'screen_classes']]]),
('3D', 
   [['Initial Model', 'rct', 'initvolume_ransac', 'preprocess_volume'], 
    ['Model Refinement', 'projmatch', 'ml3d', 'relion3d'],
    ['Analysis', ['Flexibility', 'nma', 'nma_alignment'], 'create_volume_mask']])
,
('Other',
 [['Extra', 'custom','subtraction', 'mltomo']])
]

#--------------------------------------------------------------------------------
# Protocol data support
#--------------------------------------------------------------------------------
class ProtocolData:
    def __init__(self, section, group, name, title, path):
        self.section = section
        self.group = group
        self.name = name
        self.title = title
        self.dir = path

class ProtocolDictionary(dict):
    def __init__(self):
        for section, sectionList in sections:
            for groupList in sectionList:
                group = groupList[0]
                protocolList = groupList[1:]
                for protocol in protocolList:
                    if type(protocol) != list:
                        self.addProtocol(section, group, protocol)
                    else:
                        for p in protocol[1:]:
                            self.addProtocol(section, group, p)
        # Add special 'xmipp_program'
        self.addProtocol(None, None, 'xmipp')

    def addProtocol(self, section, group, protocol):
        title, path = protocols[protocol]
        p = ProtocolData(section, group, protocol, title, path)
        setattr(self, protocol, p)
        self[protocol] = p
        

protDict = ProtocolDictionary()

#--------------------------------------------------------------------------------
# Project default settings
#--------------------------------------------------------------------------------
PROJECT_DB = '.project.sqlite'

projectDefaults = {
                   'Cfg': '.project.cfg',
                   'Db': PROJECT_DB,
                   'LogsDir': 'Logs',
                   'RunsDir': 'Runs',
                   'TmpDir': 'Tmp',
                   'RunsPrefix': 'run',
                   'TableGroups': 'groups',
                   'TableParams': 'params',
                   'TableProtocols': 'protocols',
                   'TableProtocolsGroups': 'protocols_groups',
                   'TableRuns': 'runs',
                   'TableSteps': 'steps'
                   } 

#--------------------------------------------------------------------------------
# GUI Properties
#--------------------------------------------------------------------------------

#Font
#FontName = "Helvetica"
FontName = "Verdana"
FontSize = 10

#TextColor
CitationTextColor = "dark olive green"
LabelTextColor = "black"
SectionTextColor = "blue4"

#Background Color
BgColor = "light grey"
LabelBgColor = "white"
HighlightBgColor = BgColor
ButtonBgColor = "LightBlue"
ButtonActiveBgColor = "LightSkyBlue"
EntryBgColor = "lemon chiffon" 
ExpertLabelBgColor = "light salmon"
SectionBgColor = ButtonBgColor

#Color
ListSelectColor = "DeepSkyBlue4"
BooleanSelectColor = "white"
ButtonSelectColor = "DeepSkyBlue2"

#Dimensions limits
MaxHeight = 650
MaxWidth = 800
MaxFontSize = 14
MinFontSize = 6
WrapLenght = MaxWidth - 50
