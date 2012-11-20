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
        self.addParamsLine(' [--colors <colors>]    : Colors for each plot')
        self.addParamsLine('   alias -c;')
        self.addParamsLine(' [--markers <markers>]  : Markers for each plot')
        self.addParamsLine('   alias -m;')
        self.addParamsLine(' [--style <linestyle>]  : Line style for each plot')
        self.addParamsLine('   alias -s;')
        self.addParamsLine(' [--nbins <nbin>]          : Create histogram with Y data and nbin bins')
        self.addParamsLine('   alias -n;')
        self.addParamsLine(' [--legend <location=best>] : Select where to place legend')
        self.addParamsLine('    where <location> ')
        self.addParamsLine('     best none  upper_right upper_left lower_right lower_left')
        self.addParamsLine('   alias -l;')
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
        --title "My figure" --colors "yo bo go"')
        self.addExampleLine('Plot using different line styles:', False)
        self.addExampleLine('xmipp_metadata_plot -i results.xmd -x iterationNumber  \
        -y "sigmaNoise sigmaOffset iterationNumber" \
        --title "My figure" --colors "yellow-- blue. green.."')
        self.addExampleLine('Colors, markers and style lines are described here: http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot', False)

    def getList(self, paramName, defaultValue=[None]):
        '''Return list and len of it '''
        if self.checkParam(paramName):
            paramList = self.getParam(paramName).split()
        else:
            paramList = defaultValue
        return paramList, len(paramList)
        
        
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
        
        colors, lenColors = self.getList('--colors', ['g', 'b', 'r', 'y', 'c', 'm', 'k'])
        markers, lenMarkers = self.getList('--markers')
        styles, lenStyles = self.getList('--style')
        
        if self.checkParam('--nbins'):
            nBins = self.getIntParam('--nbins')
        else:
            nBins = None
            
        
        title = self.getParam('--title')
        xplotter = XmippPlotter()
        xplotter.createSubPlot(title, xlabel, ylabel)
        
        for i, l in enumerate(ylabels):
            c = colors[i % lenColors]
            if c.startswith('('): # Convert rgb format to tuples
                #p = c[1:-1].split(',')
                #c = (int(p[0]), int(p[1]), int(p[2]))
                c = (255, 0, 0)
            m = markers[i % lenMarkers]
            if m == "none":
                m = None
            s = styles[i % lenStyles]
            if s == "none":
                s = None
            xplotter.plotMd(md, mdLabelX, str2Label(l), color=c, marker=m, linestyle=s, nbins=nBins)#if nbins is present do an histogram
        
        legendLocation = self.getParam('--legend')
        if legendLocation != 'none':
            xplotter.showLegend(ylabels, loc=legendLocation.replace('_', ' '))
        xplotter.show()

if __name__ == '__main__':
    ScriptPlotMetadata().tryRun()
