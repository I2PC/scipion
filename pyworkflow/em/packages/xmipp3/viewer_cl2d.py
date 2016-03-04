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
import pyworkflow.utils as pwutils
from pyworkflow.viewer import ProtocolViewer, WEB_DJANGO
from pyworkflow.em import *
from protocol_cl2d import XmippProtCL2D
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
from pyworkflow.protocol.params import LabelParam, BooleanParam, LEVEL_ADVANCED
import glob

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
              label='Levels selection', condition='doShowLastLevel==%d' % LEVEL_SEL,
              help='Specify a  list of levels like: 0,1,3 or 0-3 ') 

        # Sorted or not
        form.addParam('doSorted', BooleanParam, default=False, 
                      label="Display sorted classes",
                      help="Note that the average of the images assigned to the sorted class are not, in general, aligned with the sorted class, "
                           "because the sorted class representatives have been realigned to be shown all in the same orientation")     
         
        # Classes
        CLASS_CHOICES = ['classes']
        if self.protocol.doCore:
            CLASS_CHOICES.append('cores')
            if self.protocol.doStableCore:
                CLASS_CHOICES.append('stable cores')
        form.addParam('classesToShow', EnumParam, choices=CLASS_CHOICES,
                      label="What to show", default=CLASSES,
                      display=EnumParam.DISPLAY_LIST)
        # Convergence
        form.addParam('showConvergence', LabelParam,
                      label="Visualize convergence", 
                      help="It shows per iteration the similarity (cl2dSimilarity) between "
                           "experimental images and the references as well as the number "
                           "experimental images changing class.")
        form.addParam('showClassHierarchy', LabelParam,  
                      expertLevel=LEVEL_ADVANCED, hidden=not self.protocol.computeHierarchy,
                      label="Visualize class hierarchy.")      
    
    def _getVisualizeDict(self):
        return {'classesToShow': self._viewLevelFiles,
                'showConvergence': self._viewConvergence,
                'showClassHierarchy': self._viewClassHierarchy,
                }        

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

    def _viewInfo(self, fn):
        """ Display the info block in the metadata for a level convergence. """
        return DataView('info@' + fn)
    
    def _viewConvergence(self, e=None):
        levelFiles = self.protocol._getLevelMdFiles()
        views = []
        errors = []
        
        if levelFiles:
            if self.doShowLastLevel == LEVEL_LAST:
                views.append(self._viewInfo(levelFiles[-1]))
            else:
                if self.showSeveralLevels.empty():
                    errors.append('Please select the levels that you want to visualize.')
                else:
                    listOfLevels = []
                    try:
                        listOfLevels = self._getListFromRangeString(self.showSeveralLevels.get())
                    except Exception, ex:
                        errors.append('Invalid levels range.')
                        
                    files = []
                    for level in listOfLevels:
                        fn = self.protocol._getExtraPath("level_%02d/level_classes.xmd" % level)
                        if os.path.exists(fn):
                            files.append(fn)
                        else:
                            errors.append('Level %s does not exist.' % level)
                    
                    for fn in files:                        
                        views.append(self._viewInfo(fn))
        else:
            errors.append('Classes have not been produced.')
                        
        self.errorList(errors, views)
        
        return views                     
                    
    def _viewClassHierarchy(self, e=None):
        fnSubset = self._getSubset()
        fnHierarchy = self.protocol._getExtraPath("classes%s_hierarchy.txt" % fnSubset)
        if os.path.exists(fnHierarchy):
            return [self.textView([fnHierarchy])] 
            
    def _getInputParticles(self):
        return self.protocol.inputParticles.get()
    
    def _setClassOrder(self, classItem, classRow):
        """ Set the class order for display the classes sorted. """
        lastOrder = getattr(self, 'lastOrder', 0) + 1
        classItem._xmipp_cl2dOrder = Integer(lastOrder)
        setattr(self, 'lastOrder', lastOrder)
        
    def _getClassesSqlite(self, blockName, fn):
        """ Read the classes from Xmipp metadata and write as sqlite file. """
        fnSqlite = fn + '_'+blockName+'.sqlite'

        # Generate the sqlite file from classes xmd if either
        # 1) Sqlite files has not been generated or
        # 2) Xmd file is newer than sqlite (using modification time)
        if os.path.exists(fnSqlite):
            if os.path.getmtime(fn) > os.path.getmtime(fnSqlite)-10:
                pwutils.cleanPath(fnSqlite) # Clean to load from scratch

        if not os.path.exists(fnSqlite):
            classesSet = SetOfClasses2D(filename=fnSqlite)
            classesSet.setImages(self._getInputParticles())
            readSetOfClasses(classesSet, fn, blockName, 
                             preprocessClass=self._setClassOrder)
            classesSet.write()
            
        return fnSqlite
    
    def _viewClasses(self, blockName, fn, viewParams=None):
        fnSqlite = self._getClassesSqlite(blockName, fn)
        return ClassesView(self._project, self.protocol.strId(), fnSqlite,
                           other=self._getInputParticles().strId(),
                           viewParams=viewParams or {})
    
    def _viewLevelFiles(self, e=None):
        fnSubset = self._getSubset()
        views = []
        errors = []
        #obj = getattr(self.protocol, "outputClasses" + fnSubset)
        levelFiles = self.protocol._getLevelMdFiles(fnSubset)

        if levelFiles:
            suffix = ''
            viewParams = {}
            
            if self.doSorted:
                suffix = '_sorted'
                viewParams[showj.SORT_BY] = '_xmipp_cl2dOrder'
                
            if self.doShowLastLevel == LEVEL_LAST:
                views.append(self._viewClasses("classes"+suffix, levelFiles[-1], viewParams))
            else:
                if self.showSeveralLevels.empty():
                    errors.append('Please select the levels that you want to visualize.')
                else:
                    listOfLevels = []
                    try:
                        listOfLevels = self._getListFromRangeString(self.showSeveralLevels.get())
                    except Exception, ex:
                        errors.append('Invalid levels range.')
                        
                    files = []
                    for level in listOfLevels:
                        fn = self.protocol._getExtraPath("level_%02d/level_classes%s.xmd"%(level,fnSubset))
                        if os.path.exists(fn):
                            files.append(("classes"+suffix, fn))
                        else:
                            errors.append('Level %s does not exist.' % level)
                    
                    for block, fn in files:                        
                        views.append(self._viewClasses(block, fn, viewParams))
        else:
            errors.append('Classes %s have not been produced.' % fnSubset.replace("_", " "))
        
        self.errorList(errors, views)
            
        return views        

