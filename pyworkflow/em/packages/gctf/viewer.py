# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import exists

from pyworkflow.gui.project import ProjectWindow
from pyworkflow.viewer import (DESKTOP_TKINTER, WEB_DJANGO, Viewer)
import pyworkflow.em as em
import pyworkflow.em.showj as showj
import pyworkflow.utils as pwutils
from pyworkflow.em.plotter import EmPlotter
from protocol_gctf import ProtGctf
from convert import readCtfModel


OBJCMD_GCTF = "Display Ctf Analysis"


class ProtGctfViewer(Viewer):
    """ Visualization of Gctf."""
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Gctf'
    _targets = [ProtGctf]
     
    def __init__(self, **args):
        Viewer.__init__(self, **args)
        self._views = []
         
    def visualize(self, obj, **args):
        self._visualize(obj, **args)
         
        for v in self._views:
            v.show()
             
    def _visualize(self, obj, **args):
        cls = type(obj)
         
        def _getMicrographDir(mic):
            """ Return an unique dir name for results of the micrograph. """
            from pyworkflow.utils.path import removeBaseExt
            return obj._getExtraPath(removeBaseExt(mic.getFileName()))
         
        def iterMicrographs(mics):
            """ Iterate over micrographs and yield
            micrograph name and a directory to process.
            """
            for mic in mics:
                micFn = mic.getFileName()
                micDir = _getMicrographDir(mic) 
                yield (micFn, micDir, mic)
         
        def visualizeObjs(obj, setOfMics):
            if exists(obj._getPath("ctfs_temporary.sqlite")):
                os.remove(obj._getPath("ctfs_temporary.sqlite"))
             
            ctfSet = self.protocol._createSetOfCTF("_temporary")
            for fn, micDir, mic in iterMicrographs(setOfMics):
                samplingRate = mic.getSamplingRate() * self.protocol.ctfDownFactor.get()
                mic.setSamplingRate(samplingRate)
                out = self.protocol._getCtfOutPath(micDir)
                psdFile = self.protocol._getPsdPath(micDir)
                
                if exists(out) and exists(psdFile):
                    ctfModel = em.CTFModel()
                    readCtfModel(ctfModel, out)
                    ctfModel.setPsdFile(psdFile)
                    ctfModel.setMicrograph(mic)
                    ctfSet.append(ctfModel)
            
            if ctfSet.getSize() < 1:
                raise Exception("Has not been completed the CTF estimation of any micrograph")
            else:
                ctfSet.write()
                ctfSet.close()
                self._visualize(ctfSet)
        
        
        if issubclass(cls, ProtGctf) and not obj.hasAttribute("outputCTF"):
            mics = obj.inputMicrographs.get()
            visualizeObjs(obj, mics)
        elif obj.hasAttribute("outputCTF"):
            self._visualize(obj.outputCTF)
        else:
            fn = obj.getFileName()
            if obj.strId() == "None":
                objName = fn
            else:
                objName = obj.strId()
            psdLabels = '_psdFile'
            labels = 'id enabled comment %s _defocusU _defocusV _defocusAngle _defocusRatio' % psdLabels
            labels = labels + ' _gctf_ctfResolution _gctf_crossCorrelation _micObj._filename'
            print "objName, ", objName
            self._views.append(em.ObjectView(self._project, objName, fn,
                                             viewParams={showj.MODE: showj.MODE_MD,
                                                         showj.ORDER: labels,
                                                         showj.VISIBLE: labels,
                                                         showj.ZOOM: 50,
                                                         showj.RENDER: psdLabels,
                                                         showj.OBJCMDS: "'%s'" % OBJCMD_GCTF}))

        return self._views


def createCtfPlot(ctfSet, ctfId):
    ctfModel = ctfSet[ctfId]
    psdFn = ctfModel.getPsdFile()
    fn = pwutils.removeExt(psdFn) + "_EPA.txt"
    gridsize = [1, 1]
    xplotter = EmPlotter(x=gridsize[0], y=gridsize[1], windowTitle='CTF Fitting')
    plot_title = "CTF Fitting"
    a = xplotter.createSubPlot(plot_title, 'Resolution (Angstroms)', 'CTF', yformat=False)
    a.invert_xaxis()
    for i in range(1, 5):
        _plotCurve(a, i, fn)
    xplotter.showLegend(['simulated CTF',
                         'equiphase avg.',
                         'equiphase avg. - bg',
                         'cross correlation'])
    a.grid(True)
    xplotter.show()


ProjectWindow.registerObjectCommand(OBJCMD_GCTF, createCtfPlot)


def _plotCurve(a, i, fn):
    freqs = _getValues(fn, 0)
    curv = _getValues(fn, i)
    a.plot(freqs, curv)


def _getValues(fn, col):
    f = open(fn)
    values = []
    for line in f:
        if not line.startswith('Resolution', 2, 12):
            column = line.split()
            value = float(column[col])
            values.append(value)
    f.close()
    return values

