#--------------------------------------------------------------------------------
# Protocol information
#--------------------------------------------------------------------------------
class ProtocolData:
    def __init__(self,key,title,protocolDir):
        self.key=key
        self.title=title
        self.protocolDir=protocolDir

class ProtocolDictionary:
    def addProtocol(self,key,title,protocolDir):
        p = ProtocolData(key,title,protocolDir)
        self.protocolDict[key]=p
        return p
    
    def __init__(self):
        self.protocolDict={}
        self.preprocess_micrographs = self.addProtocol('preprocess_micrographs', 'Preprocess Micrograph', 'Preprocess')
        self.particle_pick = self.addProtocol('particle_pick',  'Particles picking', 'ParticlePicking')
        self.preprocess_particles = self.addProtocol('preprocess_particles',  'Preprocess Particles', 'Images')
        self.ml2d = self.addProtocol('ml2d', 'ML2D', '2D/ML2D')
        self.cl2d = self.addProtocol('cl2d', 'CL2D', '2D/CL2D')
        self.kerdensom = self.addProtocol('kerdensom', 'KerDenSOM',  '2D/KerDenSOM')
        self.rotspectra = self.addProtocol('rotspectra', 'Rotational Spectra', '2D/RotSpectra')
        self.commonlines = self.addProtocol('commonlines', 'Common Lines', '3D/CommonLines')
        self.rct = self.addProtocol('rct', 'Random Conical Tilt', '3D/RCT')
        self.projmatch = self.addProtocol('projmatch', 'Projection Matching', '3D/ProjMatch') 
        self.projsubs = self.addProtocol('subtraction', 'Partial Projection Subtraction', '3D/ProjSubs')

protDict=ProtocolDictionary()

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
                   'TableSteps': 'steps',
                   'TableStepsRestart': 'steps_restart'
                   } 

#--------------------------------------------------------------------------------
# GUI Menus
#--------------------------------------------------------------------------------

sections = [
('Preprocessing', 
   [['Preprocess Micrograph', protDict.preprocess_micrographs.key], 
    ['Particles picking', protDict.particle_pick.key], 
    ['Preprocess Particles', protDict.preprocess_particles.key]]),
('2D', 
   [['Align+Classify', protDict.ml2d.key, protDict.cl2d.key], 
    ['Align', protDict.ml2d.key, protDict.cl2d.key], 
    ['Classify', protDict.kerdensom.key, protDict.rotspectra.key]]),
('3D', 
   [['Initial Model', protDict.commonlines.key, protDict.rct.key], 
    ['Model Refinement', protDict.projmatch.key]])
,
('Other',
 [['Extra',protDict.projsubs.key]])
]

def getSectionByKey(prot): 
    for s, list in sections:
        ss = []
        for subList in list:
            if prot.key in subList:
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
BgColor = "white"
LabelBgColor = BgColor
HighlightBgColor = BgColor
ButtonBgColor = "LightBlue"
ButtonActiveBgColor = "LightSkyBlue"
EntryBgColor = "lemon chiffon" 
ExpertLabelBgColor = "light salmon"

#Color
ListSelectColor = "DeepSkyBlue4"
BooleanSelectColor = "white"

#Dimensions limits
MaxHeight = 800
MaxWidth = 800
MaxFontSize = 14
MinFontSize = 6
