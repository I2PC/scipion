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
        protocolDict[key]=p
        return p
    protocolDict={}
    preprocess_micrographs = addProtocol('preprocess_micrographs', 'Preprocess Micrograph', 'Preprocess')
    particle_pick = addProtocol('particle_pick',  'Particles picking', 'ParticlePicking')
    preprocess_particles = addProtocol('preprocess_particles',  'Preprocess Particles', 'Images')
    ml2d = addProtocol('ml2d', 'ML2D', '2D/ML2D')
    cl2d = addProtocol('cl2d', 'CL2D', '2D/CL2D')
    kerdensom = addProtocol('kerdensom', 'KerDenSOM',  '2D/KerDenSOM')
    rotspectra = addProtocol('rotspectra', 'Rotational Spectra', '2D/RotSpectra')
    commonlines = addProtocol('commonlines', 'Common Lines', '3D/CommonLines')
    rct = addProtocol('rct', 'Random Conical Tilt', '3D/RCT')
    projmatch = addProtocol('projmatch', 'Projection Matching', '3D/ProjMatch') 
    projsubs = addProtocol('subtraction', 'Partial Projection Subtraction', '3D/ProjSubs')

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
   [['Preprocess Micrograph', ProtocolDictionary.preprocess_micrographs], 
    ['Particles picking', ProtocolDictionary.particle_pick], 
    ['Preprocess Particles', ProtocolDictionary.preprocess_particles]]),
('2D', 
   [['Align+Classify', ProtocolDictionary.ml2d, ProtocolDictionary.cl2d], 
    ['Align', ProtocolDictionary.ml2d, ProtocolDictionary.cl2d], 
    ['Classify', ProtocolDictionary.kerdensom, ProtocolDictionary.rotspectra]]),
('3D', 
   [['Initial Model', ProtocolDictionary.commonlines, ProtocolDictionary.rct], 
    ['Model Refinement', ProtocolDictionary.projmatch]]),
('Other',
 [['Browse',ProtocolDictionary.projsubs]])
]

def getSectionByKey(protKey):
    for s, list in sections:
        for subList in list:
            if protKey in subList:
                return (s, subList[0])
    return None

def getSectionByValue(protValue):
    for k, v in ProtocolDictionary.protocolDict.iteritems():
        if v == protValue:
            return getSectionByKey(k)

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
BooleanSelectColor = "DeepSkyBlue4"

#Dimensions limits
MaxHeight = 800
MaxWidth = 800
MaxFontSize = 14
MinFontSize = 6
