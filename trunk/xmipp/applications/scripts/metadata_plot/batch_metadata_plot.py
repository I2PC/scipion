#!/usr/bin/env xmipp_python
"""/***************************************************************************
 *
 * Authors:     Roberto Marabini
 *              J. M. de la Rosa Trevin
 *
 * Universidad Autonoma de Madrid
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
"""
#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import XmippScript

class ScriptPlotMetadata(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Plot some values from metadata.')
        ## params
        self.addParamsLine(' -i <metadata>          : Input metadata')
        self.addParamsLine('   alias --input;')
        self.addParamsLine(' [-x <label>]           : Label for X axis data')
        self.addParamsLine('   alias --xlabel;')
        self.addParamsLine('-y <labels>             : Labels list for Y data')
        self.addParamsLine('   alias --ylabel;')
        self.addParamsLine(' [--title <title="">]   : Plot title')
        self.addParamsLine(' [--xtitle <title>]     : Plot x axis label')
        self.addParamsLine(' [--ytitle <title>]     : Plot y axis label')
        self.addParamsLine(' [-c <colors>]          : Colors for each plot')
        self.addParamsLine('   alias --colors;')
        self.addParamsLine(' [-n <nbin>]          : Create histogram with Y data and nbin bins')
        self.addParamsLine('   alias --nbins;')
        ## examples
        self.addExampleLine('Simple plot of label "sigmaNoise" from metadata', False)
        self.addExampleLine('xmipp_metadata_plot -i results.xmd -y sigmaNoise')
        self.addExampleLine('Additionally take values for X from other label and set title', False)
        self.addExampleLine('xmipp_metadata_plot -i results.xmd -x iterationNumber -y sigmaNoise --title "My figure" ')
        self.addExampleLine('Plot different labels and select colors:', False)
        self.addExampleLine('xmipp_metadata_plot -i results.xmd -x iterationNumber -y "sigmaNoise sigmaOffset iterationNumber" --title "My figure" --colors "yellow blue green"')
        self.addExampleLine('Plot using dots:', False)
        self.addExampleLine('xmipp_metadata_plot -i results.xmd -x iterationNumber  \
        -y "sigmaNoise sigmaOffset iterationNumber" \
        --title "My figure" --colors "yellowo blueo greeno"')
        self.addExampleLine('Plot using different line styles:', False)
        self.addExampleLine('xmipp_metadata_plot -i results.xmd -x iterationNumber  \
        -y "sigmaNoise sigmaOffset iterationNumber" \
        --title "My figure" --colors "yellow-- blue. green.."')
        self.addExampleLine('Colors, markers and style lines are described here: http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot', False)

    def run(self):        
        from xmipp import MetaData, str2Label
        from protlib_gui_figure import XmippPlotter        
        
        md = MetaData(self.getParam('-i'))
        if self.checkParam('--xlabel'):
            xlabel = self.getParam('--xlabel')
            mdLabelX = str2Label(xlabel)
        else:
            xlabel = ""
            mdLabelX = None
            
        if self.checkParam('--xtitle'):
            xlabel = self.getParam('--xtitle')
        ylabels = self.getParam('--ylabel').split()
        if self.checkParam('--ytitle'):
            ylabel = self.getParam('--ytitle')
        else:
            ylabel = ylabels[0]
        if self.checkParam('--colors'):
            colors = self.getParam('--colors').split()
        else:
            colors = ['g', 'b', 'r', 'y']    
        if self.checkParam('--nbins'):
            nBins = int(self.getParam('--nbins'))
        else:
            nBins = None
            
        
        title = self.getParam('--title')
        xplotter = XmippPlotter()
        xplotter.createSubPlot(title, xlabel, ylabel)
        
        for i, l in enumerate(ylabels):
            if nBins:
                xplotter.plotMd(md, mdLabelX, str2Label(l), color=colors[i], nbins=nBins)#if nbins presnts do an histogram
            else:
                xplotter.plotMd(md, mdLabelX, str2Label(l), color=colors[i])#if nbins presnts do an histogram
        xplotter.showLegend(ylabels)
        xplotter.show()

if __name__ == '__main__':
    ScriptPlotMetadata().tryRun()
