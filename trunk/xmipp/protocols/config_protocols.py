#--------------------------------------------------------------------------------
# Protocol information
#--------------------------------------------------------------------------------
class ProtocolData:
    def __init__(self,name,title,dir):
        self.name=name
        self.title=title
        self.dir=dir

class ProtocolDictionary(dict):
    def addProtocol(self,name,title,dir):
        p = ProtocolData(name,title,dir)
        self[name]=p
        return p
    
    def __init__(self):
        self.import_micrographs = self.addProtocol('import_micrographs', 'Import Micrographs', 'ImportedMicrographs')
        self.particle_pick = self.addProtocol('particle_pick',  'Manual picking', 'ParticlePicking/Manual')
        self.particle_pick_auto = self.addProtocol('particle_pick_auto',  'Automatic picking', 'ParticlePicking/Auto')
        self.extract_particles = self.addProtocol('extract_particles',  'Extract Particles', 'Images/Extracted')
        self.import_particles = self.addProtocol('import_particles',  'Import Particles', 'Images/Imported')
        self.preprocess_particles = self.addProtocol('preprocess_particles',  'Preprocess Particles', 'Images/Preprocessed')
        self.ml2d = self.addProtocol('ml2d', 'ML2D', '2D/ML2D')
        self.cl2d = self.addProtocol('cl2d', 'CL2D', '2D/CL2D')
        self.cl2d_alignment = self.addProtocol('cl2d_align', 'Only align', '2D/Alignment')
        self.kerdensom = self.addProtocol('kerdensom', 'KerDenSOM',  '2D/KerDenSOM')
        self.rotspectra = self.addProtocol('rotspectra', 'Rotational Spectra', '2D/RotSpectra')
        self.commonlines = self.addProtocol('commonlines', 'Common Lines', '3D/CommonLines')
        self.rct = self.addProtocol('rct', 'Random Conical Tilt', '3D/RCT')
        self.projmatch = self.addProtocol('projmatch', 'Projection Matching', '3D/ProjMatch') 
        self.ml3d = self.addProtocol('ml3d', 'ML3D', '3D/ML3D')
        self.projsubs = self.addProtocol('subtraction', 'Partial Projection Subtraction', '3D/ProjSubs')
        self.dummy = self.addProtocol('dummy', 'Dummy', 'Dummy')
        self.xmipp_program = self.addProtocol('xmipp', 'Xmipp Programs', 'XmippPrograms')
protDict = ProtocolDictionary()

#--------------------------------------------------------------------------------
# Project default settings
#--------------------------------------------------------------------------------
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
                   'TableSteps': 'steps'
                   } 

#--------------------------------------------------------------------------------
# GUI Menus
#--------------------------------------------------------------------------------

sections = [
('Preprocessing', 
   [['Import Micrographs', protDict.import_micrographs.name], 
    ['Particle picking', protDict.particle_pick.name, protDict.particle_pick_auto.name, protDict.extract_particles.name, protDict.import_particles.name], 
    ['Preprocess Particles', protDict.preprocess_particles.name]]),
('2D', 
   [['Align+Classify', protDict.ml2d.name, protDict.cl2d.name, protDict.cl2d_alignment.name], 
    ['Classify', protDict.kerdensom.name, protDict.rotspectra.name]]),
('3D', 
   [['Initial Model', protDict.commonlines.name, protDict.rct.name], 
    ['Model Refinement', protDict.projmatch.name, protDict.ml3d.name]])
,
('Other',
 [['Extra',protDict.projsubs.name, protDict.dummy.name]])
]

def getSectionByName(prot): 
    for s, list in sections:
        ss = []
        for subList in list:
            if prot.name in subList:
                ss.append(subList[0])
        if len(ss) > 0:
            return (s, ss)
    return None


#--------------------------------------------------------------------------------
# GUI Properties
#--------------------------------------------------------------------------------

#Font
FontName = "Helvetica"
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
ButtonSelectColor = "DodgerBlue3"

#Dimensions limits
MaxHeight = 650
MaxWidth = 800
MaxFontSize = 14
MinFontSize = 6
WrapLenght = MaxWidth - 50