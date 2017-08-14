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

from pyworkflow.gui.project import ProjectWindow
import pyworkflow.utils as pwutils
from pyworkflow.em.plotter import EmPlotter



def createCtfPlot(ctfSet, ctfId):
    ctfModel = ctfSet[ctfId]
    psdFn = ctfModel.getPsdFile()
    fn = pwutils.removeExt(psdFn) + "_EPA.txt"
    gridsize = [1, 1]
    xplotter = EmPlotter(x=gridsize[0], y=gridsize[1],
                         windowTitle='CTF Fitting')
    plot_title = "CTF Fitting"
    a = xplotter.createSubPlot(plot_title, 'Resolution (Angstroms)', 'CTF',
                               yformat=False)
    a.invert_xaxis()
    for i in range(1, 5):
        _plotCurve(a, i, fn)
    xplotter.showLegend(['simulated CTF',
                         'equiphase avg.',
                         'equiphase avg. - bg',
                         'cross correlation'])
    a.grid(True)
    xplotter.show()


OBJCMD_GCTF = "Display Ctf Analysis"

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

