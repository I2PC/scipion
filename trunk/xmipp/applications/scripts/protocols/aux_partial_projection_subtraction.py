
from xmipp import *
def joinImageCTF(_log, dict):
    tmpMD = MetaData(dict['CTFgroupName'])
    MDaux = MetaData(dict['filename_currentAngles'])
    outMD = MetaData()
    outMD.join(MDaux, tmpMD, MDL_IMAGE, NATURAL_JOIN)
    outMD.write(dict['DocFileRef'], MD_APPEND)

