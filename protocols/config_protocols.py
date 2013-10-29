
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
        'kerdensom': ('KerDenSOM', '2D/KerDenSOM'),
        'rotspectra': ('Rotational Spectra', '2D/RotSpectra'),
        'screen_classes': ('Screen classes', '2D/Screening'),
        'rct': ('Random Conical Tilt', '3D/InitialVolume/RCT'),
        'initvolume_ransac': ('RANSAC', '3D/InitialVolume/RANSAC'),
        'initvolume_simanneal': ('Simulated Annealing', '3D/InitialVolume/SimAnneal'),
        'initvolume_validation': ('Validate Volume', '3D/InitialVolume/InitVolumeValidation'),
        'convert_pdb': ('Convert PDB', '3D/PDB'),
        'preprocess_volume': ('Preprocess', '3D/Preprocessed'),
        'create_volume_mask': ('Create mask', '3D/Mask'),
        'projmatch': ('Projection Matching', '3D/ProjMatch'),
        'ml3d': ('ML3D', '3D/ML3D'),
        'hg3d': ('HG3D', '3D/InitialVolume/HG3D'),
        'nma': ('Normal Mode Analysis', '3D/NMA'),
        'nma_alignment': ('Flexible alignment', '3D/NMA_alignment'),
        'resolution3D': ('Resolution 3D', '3D/Resolution'),
        'align_volume': ('Align Volume', '3D/AlignVolume'),
        'helical_params': ('Helical Parameters', '3D/Helical'),
        'relion_classify': ('3D Classification ', '3D/RelionClass'),
        'relion_refine': ('Angle Refinement ', '3D/RelionRef'),
        'cltomo': ('CLTomo', '3D/CLTomo'),
        'mltomo': ('MLTomo', '3D/MLTomo'),
        'subtraction': ('Partial Projection Subtraction', '3D/ProjSub'),
        'custom': ('Custom', 'Custom'),
        'image_operate': ('Image Operate', 'Tools/ImageOperate'),
        'metadata_utilities': ('Metadata Utilities', 'Tools/MetadataUtilities'),
        'metadata_split': ('Metadata Split', 'Tools/MetadataSplit'),
        # 'xmipp': ('Xmipp Programs', 'XmippPrograms'), 
        'emx_import_micrographs': ('Import micrographs', 'Micrographs/EmxImported'),
        'emx_import_particles': ('Import particles', 'Images/EmxImported'),
        'emx_export_micrographs': ('Export micrographs', 'Micrographs/EmxExported'),
        'emx_export_particles': ('Export particles', 'Images/EmxExported'),
        'xray_import': ('Import tomograms', 'XrayTomo/Imported'),
        'xray_fast_align': ('Fast align tomograms', 'XrayTomo/FastAlignment'),
        'xray_reconstruct': ('Reconstruct tomograms', 'XrayTomo/Reconstruct')
        }

#--------------------------------------------------------------------------------
# Protocols sections and groups
#--------------------------------------------------------------------------------
sections = [
('Preprocessing',
   [['Micrographs', 'import_micrographs', 'screen_micrographs', 'downsample_micrographs'],
    ['Particle picking', 'particle_pick', 'particle_pick_auto'],
    ['Particles', 'extract_particles', 'import_particles', ['Other', 'preprocess_particles', 'screen_particles', 'merge_particles']]]),
('2D',
   [['Align+Classify', 'cl2d', 'ml2d', ['Other', 'cl2d_align', 'kerdensom', 'rotspectra', 'screen_classes']]]),
('3D',
   [['Initial Model', 'rct', 'initvolume_ransac', 'initvolume_simanneal', 'convert_pdb',['Heterogeneity', 'hg3d'],'initvolume_validation'],
    ['Model Refinement', 'projmatch', 'ml3d', ['relion', 'relion_classify', 'relion_refine']],
    ['Volumes', 'create_volume_mask', 'preprocess_volume', 'resolution3D', 'align_volume', 'helical_params']]),
('Other',
 [['Extra',
   'custom',
   ['Flexibility', 'nma', 'nma_alignment'],
   ['Virus', 'subtraction'],
   ['Tomography', 'mltomo', 'cltomo'],
   ['X-ray', 'xray_import', 'xray_fast_align', 'xray_reconstruct'],
   ['Tools', 'image_operate', 'metadata_utilities', 'metadata_split'],
   ['EMX', 'emx_import_micrographs', 'emx_import_particles', 'emx_export_micrographs', 'emx_export_particles']]])
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
        self.protocolPaths = []
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
        # self.addProtocol(None, None, 'xmipp')
    def existsPrefix(self, path):
        """ Find if another protocol contains this path as prefix"""
        path += '/'
        for p in self.protocolPaths:
            if p.startswith(path):
                return True
        return False

    def addProtocol(self, section, group, protocol):
        title, path = protocols[protocol]
        if self.existsPrefix(path):
            raise Exception("Path is already existing as prefix: " + path)
        self.protocolPaths.append(path)
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

# Font
# FontName = "Helvetica"
# Try to read FontName and FontSize 
# from environment variables
import os
FontName = os.environ.get('XMIPP_FONT_NAME', "Verdana")
FontSize = int(os.environ.get('XMIPP_FONT_SIZE', 10))

# TextColor
CitationTextColor = "dark olive green"
LabelTextColor = "black"
SectionTextColor = "blue4"

# Background Color
BgColor = "light grey"
LabelBgColor = "white"
HighlightBgColor = BgColor
ButtonBgColor = "LightBlue"
ButtonActiveBgColor = "LightSkyBlue"
EntryBgColor = "lemon chiffon" 
ExpertLabelBgColor = "light salmon"
SectionBgColor = ButtonBgColor

# Color
ListSelectColor = "DeepSkyBlue4"
BooleanSelectColor = "white"
ButtonSelectColor = "DeepSkyBlue2"

# Dimensions limits
MaxHeight = 650
MaxWidth = 2048
MaxFontSize = 18
MinFontSize = 10
WrapLenght = MaxWidth - 50
