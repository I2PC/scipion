# **************************************************************************
# *
# * Authors:     L. del Cano (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pyworkflow.em.packages.xmipp3.convert import readSetOfClasses
"""
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""

import os
import pyworkflow.em as em
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.protocol.params import LabelParam, BooleanParam, LEVEL_ADVANCED
from protocol_cl2d import XmippProtCL2D
import pyworkflow.em.showj as showj
from pyworkflow.protocol.params import EnumParam, StringParam


CLASSES = 0
CLASS_CORES = 1
CLASS_STABLE_CORES = 2   

LEVEL_LAST = 0
LEVEL_SEL = 1
LEVEL_CHOICES = ['last', 'selection']
        
        
class XmippCL2DViewer(ProtocolViewer):
    """ Visualization of Xmipp CL2D protocol results. """

    _label = 'viewer cl2d'
    _targets = [XmippProtCL2D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        # Select the level to show
        form.addParam('doShowLastLevel', EnumParam, default=LEVEL_LAST, 
                      choices=LEVEL_CHOICES,
                      display=EnumParam.DISPLAY_HLIST,
                      label="Level to visualize")     
        form.addParam('showSeveralLevels', StringParam, default='',
                      label='Levels selection',
                      condition='doShowLastLevel==%d' % LEVEL_SEL,
                      help='Specify a  list of levels like: 0,1,3 or 0-3 ')
        # TO-DO; All commented code is for revision with @cossorzano
        
        # Sorted or not
        # form.addParam('doSorted', BooleanParam, default=False,
        #               label="Display sorted classes",
        #               help="Note that the average of the images assigned to "
        #                    "the sorted class are not, in general, aligned with "
        #                    "the sorted class, because the sorted class "
        #                    "representatives have been realigned to be shown "
        #                    "all in the same orientation")
        # Classes
        CLASS_CHOICES = ['classes']
        if self.protocol.doCore:
            CLASS_CHOICES.append('cores')
            if self.protocol.doStableCore:
                CLASS_CHOICES.append('stable cores')
        form.addParam('classesToShow', EnumParam, choices=CLASS_CHOICES,
                      label="What to show", default=CLASSES,
                      display=EnumParam.DISPLAY_HLIST)
        # Convergence
        form.addParam('showConvergence', LabelParam,
                      label="Visualize convergence", 
                      help="It shows per iteration the similarity "
                           "(cl2dSimilarity) between experimental images and "
                           "the references as well as the number experimental "
                           "images changing class.")
        form.addParam('showClassHierarchy', LabelParam,  
                      expertLevel=LEVEL_ADVANCED,
                      hidden=not self.protocol.computeHierarchy,
                      label="Visualize class hierarchy.")

    def _getVisualizeDict(self):
        self.protocol._defineFileNames()  # Load filename templates
        return {'classesToShow': self._viewLevelFiles,
                'showConvergence': self._viewConvergence,
                'showClassHierarchy': self._viewClassHierarchy,
                }

# ==============================================================================
# showClasses
# ==============================================================================
    def _viewLevelFiles(self, e=None):
        views = []
        fnSubset = self._getSubset()
        viewFinalClasses = False
        
        # TO-DO; All commented code is for revision with @cossorzano
        viewParams = {}
    
        # if self.doSorted:
        #     suffix = '_sorted'
        #     viewParams[showj.SORT_BY] = '_xmipp_cl2dOrder'
        
        if self.doShowLastLevel == LEVEL_LAST:
            fn = self.protocol._getFileName("final_classes", sub=fnSubset)
            if os.path.exists(fn):
                print "Esto es cierto, ", fn
                viewFinalClasses = True
            levList = [self.protocol._lastLevel()]
        else:
            levList = self._getListFromRangeString(self.showSeveralLevels.get())

        for level in levList:
            if viewFinalClasses:
                fn = self.protocol._getFileName("final_classes", sub=fnSubset)
            else:
                fn = self.protocol._getFileName("level_classes", level=level,
                                                                 sub=fnSubset)
                if level <= self.protocol._lastLevel() and os.path.exists(fn):
                    fn = self.protocol._getLevelClasses(level, fnSubset)
            
            if os.path.exists(fn):
                inputParts = self.protocol.inputParticles.get()
                views.append(em.ClassesView(self.getProject(),
                                            self.protocol.strId(), fn,
                                            inputParts.strId()))
        return  views

# ==============================================================================
# showConvergence
# ==============================================================================
    def _viewConvergence(self, e=None):
        views = []
        if self.doShowLastLevel == LEVEL_LAST:
            fn = self.protocol._getFileName('level_classes',
                                            level=self.protocol._lastLevel(),
                                            sub='')
            views.append(self._viewInfo(fn))
        else:
            levList = self._getListFromRangeString(self.showSeveralLevels.get())
            if len(levList) > 0:
                for level in levList:
                    if level <= self.protocol._lastLevel():
                        fn = self.protocol._getFileName('level_classes',
                                                        level=level, sub='')
                        views.append(self._viewInfo(fn))
        return views

    def _viewInfo(self, fn):
        """ Display the info block in the metadata for a level convergence. """
        return em.DataView('info@' + fn)
        
# ==============================================================================
# showClassHierarchy
# ==============================================================================
    def _viewClassHierarchy(self, e=None):
        fnSubset = self._getSubset()
        fnHierarchy = self.protocol._getFileName('classes_hierarchy',
                                                 sub=fnSubset)
        if os.path.exists(fnHierarchy):
            return [self.textView([fnHierarchy])]

#===============================================================================
# Utils Functions
#===============================================================================
    def _getSubset(self):
        fnSubset = ""
        classesToShow = self.classesToShow.get()
        if classesToShow == CLASSES:
            fnSubset = ""
        elif classesToShow == CLASS_CORES:
            fnSubset = "_core"
        elif classesToShow == CLASS_STABLE_CORES:
            fnSubset = "_stable_core"
        return fnSubset

    # TO-DO; All commented code is for revision with @cossorzano
    
    # def _setClassOrder(self, classItem, classRow):
    #     """ Set the class order for display the classes sorted. """
    #     lastOrder = getattr(self, 'lastOrder', 0) + 1
    #     classItem._xmipp_cl2dOrder = Integer(lastOrder)
    #     setattr(self, 'lastOrder', lastOrder)

