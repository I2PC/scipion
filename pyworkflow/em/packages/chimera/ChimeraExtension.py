"""load the script in chimera for example chimera --script ScipioChimeraExt/ChimeraExtension.py
from chimera command line type one of the following three options
scipionwrite
scipionwrite model #n
scipionwrite model #n refmodel #p
"""
def cmd_scipionWrite(lowpass,args):
  from Midas.midas_text import doExtensionFunc

  def scipionWrite(model="#1",refmodel="#0"):
     #get model (pdb) id
     modelId = int(model[1:])# model to write
     refModelId = int(refmodel[1:])# coordenates system refers to this model
     
     # get actual models
     model    = chimera.openModels.list()[modelId]
     refModel = chimera.openModels.list()[refModelId]
     
     #write the model
     #wrapper should incorporate this pdb as scipion object
     from Midas import write
     write(model, refModel, "extra/model.pdb")

  doExtensionFunc(scipionWrite,args)

from Midas.midas_text import addCommand
addCommand('scipionwrite', cmd_scipionWrite, help="http://scipion.cnb.csic.es")
#TODO help
