#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------

# Author: Carlos Oscar, August 2013
#
from protlib_base import *
from os.path import exists
from protlib_utils import runJob

class ProtStructureFactor(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.structure_factor.name, scriptname, project)
        self.Import = 'from protocol_structure_factor import *'
        self.fnStructure = self.workingDirPath('structureFactor.xmd')

    def defineSteps(self):
        self.insertStep("calculateStructureFactor",Structure=self.fnStructure,InModel=self.InModel,
                        Sampling=self.Ts)

    def summary(self):
        messages = []      
        messages.append("Input: [%s]" % self.InModel)
        messages.append("Output: [%s]" % self.fnStructure)
        return messages

    def visualize(self):
        if os.path.exists(self.fnStructure):
            if self.DisplayStructureFactor:
                os.system('xmipp_metadata_plot -i %s -x resolutionFreqFourier -y resolutionLogStructure --title "Structure factor" --xtitle "Frequency (1/A)" --ytitle "Log(StructureFactor)" &'%self.fnStructure)
            if self.DisplayGuinier:
                if self.UseMatlab:
                    os.system("matlab -r \"xmipp_show_structure_factor(\'"+self.WorkingDir+"\')\"")
                else:
                    os.system('xmipp_metadata_plot -i %s -x resolutionFreqFourier2 -y resolutionLogStructure --title "Guinier plot" --xtitle "Frequency (1/A^2)" --ytitle "Log(StructureFactor)" &'%self.fnStructure)

def calculateStructureFactor(log,Structure,InModel,Sampling):
    runJob(log,"xmipp_volume_structure_factor","-i %s -o %s --sampling %f"%(InModel,Structure,float(Sampling)))
