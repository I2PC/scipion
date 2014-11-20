# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import sys
import os
from os.path import join

import Tkinter as tk
import pyworkflow.gui as gui
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider, HotButton

from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.wizard import EmWizard
from protocol_resmap import ProtResMap
from pyworkflow.gui.matplotlib_image import FigureFrame

#===============================================================================
# Pre-whitening for Resmap
#===============================================================================

class ResmapPrewhitenWizard(EmWizard):
    _targets = [(ProtResMap, ['prewhitenAng'])]

    def show(self, form):
        self.prot = form.protocol
        self.form = form
        
        emptyInput = ((not self.prot.inputVolume.get().hasValue()) or
                      (self.prot.useSplitVolume == True and 
                       not self.prot.splitVolume.get().hasValue()))
        
        if emptyInput:
            dialog.showWarning("Empty input", "Select input volume(s) first", self.form.root)
        else:
            d = PreWhiteningDialog(self.form, os.getcwd())
            if d.resultYes():
                print "result: ", d.getElbowValue(), d.getRampValue()
                form.setVar('prewhitenAng', d.getElbowValue())
                form.setVar('prewhitenRamp', d.getRampValue())


# Change default instructions message

INSTRUCTIONS = """Please check that the green line
is as straight as possible,
at least in the high frequencies.

If not, adjust the sliders below
and press the Update button.

ResMap will try to pre-whiten 
the volume again.

If you are satisfied please
press OK to use that values.    
    """    

class PreWhiteningDialog(dialog.Dialog):
    def __init__(self, form, workingDir, **kwargs):
        """ 
        Params:
            parent: parent windows of the dialog.
            provider: the TreeProvider to populate items tree.
        """
        self.form = form
        self.workingDir = workingDir
        buttons = [('Select', dialog.RESULT_YES), ('Cancel', dialog.RESULT_CANCEL)]
        dialog.Dialog.__init__(self, form.root, "Pre-whitening", buttons=buttons, default='Select', **kwargs)

    def _runBeforePreWhitening(self):
        prot = self.form.protocol
        # Convert input volumes
        ih = ImageHandler()
        ih.convert(prot.inputVolume.get(), join(self.workingDir, 'volume1.map'))
        if prot.useSplitVolume:
            ih.convert(prot.splitVolume.get(), join(self.workingDir, 'volume2.map'))
    
        self.results = prot.runResmap(self.workingDir, wizardMode=True)
        
    def _runPreWhitening(self, newElbowAngstrom, newRampWeight):
        # Add resmap libraries to the path
        sys.path.append(os.environ['RESMAP_HOME'])
        from ResMap_spectrumTools import preWhitenCube, createPrewhiteningFigure, preWhitenVolumeSoftBG
        results = self.results
        
        n = results['n']
        subVolLPF = results['subVolLPF']
        splitVolume = results['splitVolume']
        vxSize = results['vxSize']
        
        # Here we need to duplicate some code because
        # have diffent path (if/else) in the original ResMap code
        if n > subVolLPF:
            dataSize = results['cubeSize']
            preWhiteningResult = preWhitenCube( n = dataSize,
                                        vxSize        = vxSize,
                                        elbowAngstrom = newElbowAngstrom,
                                        rampWeight    = newRampWeight,
                                        dataF         = results['dataF'],
                                        dataBGF       = results['dataBGF'],
                                        dataBGSpect   = results['dataBGSpect'])
        else:
            dataSize = n
            if splitVolume == False:
                preWhiteningResult = preWhitenVolumeSoftBG( n = n,
                                        vxSize        = vxSize,
                                        elbowAngstrom = newElbowAngstrom,
                                        rampWeight    = newRampWeight,
                                        dataF         = results['dataF'],
                                        dataBGSpect   = results['dataBGSpect'],
                                        softBGmask    = results['softBGmask'],
                                        )
            else:
                preWhiteningResult = preWhitenCube( n = n,
                                        vxSize        = vxSize,
                                        elbowAngstrom = newElbowAngstrom,
                                        rampWeight    = newRampWeight,
                                        dataF         = results['dataF'],
                                        dataBGF       = results['dataBGF'],
                                        dataBGSpect   = results['dataBGSpect'])

        dataPW = preWhiteningResult['dataPW']          
        
        self.figure.clear()
        createPrewhiteningFigure(figure = self.figure,
                                 showButtons = False,
                                 showSliders = False,
                                 instructions = INSTRUCTIONS,
                                 elbowAngstrom = newElbowAngstrom,
                                 rampWeight    = newRampWeight,
                                 dataSpect     = results['dataPowerSpectrum'],
                                 dataBGSpect   = results['dataBGSpect'],
                                 peval         = preWhiteningResult['peval'],
                                 dataPWSpect   = preWhiteningResult['dataPWSpect'],
                                 dataPWBGSpect = preWhiteningResult['dataPWBGSpect'],
                                 vxSize        = vxSize,
                                 dataSlice     = results['data'][int(dataSize/2),:,:],
                                 dataPWSlice   = dataPW[int(dataSize/2),:,:]
                                )
        
    def body(self, bodyFrame):
        figFrame = FigureFrame(bodyFrame, figsize=(18, 9))
        figFrame.grid(row=0, column=0, columnspan=5)
        
        #self._runBeforePreWhitening(self.prot)
        dialog.FlashMessage(self.form.root, "Running Pre-Whitening tool...", 
                            func=self._runBeforePreWhitening)
        results = self.results
        
        self.figure = figFrame.getFigure() #plt.figure(figsize=(18, 9))
        self._runPreWhitening(results['newElbowAngstrom'], 
                              results['newRampWeight'])
        
        #bodyFrame.config()
        bodyFrame.columnconfigure(0, weight=1)
        bodyFrame.rowconfigure(0, weight=1)
        
        controlsFrame = tk.Frame(bodyFrame)
        controlsFrame.grid(row=1, column=0)
        
        self.elbowSlider = LabelSlider(controlsFrame, "Angstroms", from_=2.1*results['vxSize'], to=100, 
                             value=results['newElbowAngstrom'])
        self.elbowSlider.grid(row=1, column=0, padx=5, pady=5)
        
        self.rampSlider = LabelSlider(controlsFrame, "Ramp weight", from_=0.0, to=1., 
                             value=results['newRampWeight'])
        self.rampSlider.grid(row=1, column=1, padx=5, pady=5)
        
        self.updateBtn = HotButton(controlsFrame, "   Update   ", 
                                   command=self._onUpdate,
                                   tooltip="Update plots with new pre-whitening parameters.")
        self.updateBtn.grid(row=1, column=2, padx=10, pady=5)
        
    def getElbowValue(self):
        return self.elbowSlider.get()
        
    def getRampValue(self):
        return self.rampSlider.get()
    
    def _onUpdate(self, e=None):
        self._runPreWhitening(self.getElbowValue(), self.getRampValue())
