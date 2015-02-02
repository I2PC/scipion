# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.object import Float, Object, String, Float, Integer
from pyworkflow.utils.path import cleanPath
from pyworkflow.protocol.params import PointerParam

from pyworkflow.em import EMSet
from protocol_2d import ProtAlign2D

# For the viewer part
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView


class ProtClassesConsensus(ProtAlign2D):
    """ Compare two SetOfClasses 
    """
    _label = 'classes consensus'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputClasses1', PointerParam, pointerClass='SetOfClasses',
                      label='Input classes 1',
                      help='')
        form.addParam('inputClasses2', PointerParam, pointerClass='SetOfClasses',
                      label='Input classes 2',
                      help='')

        form.addParallelSection(threads=0, mpi=0)

#--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self._insertFunctionStep('compareClassesStep', 
                                 self.inputClasses1.getObjId(),
                                 self.inputClasses2.getObjId())

    def compareClassesStep(self, i1, i2):
        set1 = self.inputClasses1.get()
        set2 = self.inputClasses2.get()
        
        # Compare each pair of class from set1 and set2
        # compute the Jaccard index for each (J = len(intersection) / len(union))
        # Create a list will all pairs indexes and the sort them
        jaccardList = []
        f = open(self._getPath('jaccard.txt'), 'w')
        f.write('; class1 class2 intersection(i) union(i) jaccard index = len(i)/len(u)\n')
        for cls1 in set1:
            ids1 = cls1.getIdSet()
            for cls2 in set2:
                ids2 = cls2.getIdSet()
                inter = len(ids1.intersection(ids2))
                union = len(ids1.union(ids2))
                jaccardIndex = float(inter) / union
                jaccardTuple = (cls1.getObjId(), cls2.getObjId(), inter, union, jaccardIndex)
                f.write('%d %d %d %d %0.3f\n' % jaccardTuple)
                jaccardList.append(jaccardTuple)
        f.close()

        jaccardList.sort(key=lambda e: e[4], reverse=True)
        visitedClasses = set()
        outputFn = self._getPath('consensus.sqlite')
        cleanPath(outputFn)
        outputSet = EMSet(filename=outputFn)
        
        for clsId1, clsId2, inter, union, jaccardIndex in jaccardList:
            if clsId1 not in visitedClasses:
                visitedClasses.add(clsId1) # mark as visited
                cls1 = set1[clsId1]
                cls2 = set2[clsId2]
                o = Object()
                o.setObjLabel('classes %d - %d' % (clsId1, clsId2))
                o.class1 = cls1.clone()
                o.class1.id = Integer(clsId1)
                o.class2 = cls2.clone()
                o.class2.id = Integer(clsId2)
                o.jaccard = Float(jaccardIndex)
                o.intersection = Integer(inter)
                o.union = Integer(union)
                outputSet.append(o)
                
        self._defineOutputs(outputConsensus=outputSet)

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = [ ]
        return errors
    
    
class ViewerClassesConsensus(Viewer):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtClassesConsensus]
    
    def _visualize(self, obj, **kwargs):
        labels = 'class1.id class1._representative._filename class2.id class2._representative._filename jaccard intersection union'
        return [DataView(obj.outputConsensus.getFileName(), 
                         viewParams={'order': labels, 'mode': 'metadata',
                                     'visible': labels,
                                     'render': 'class1._representative._filename class2._representative._filename'
                                     })
                ]
    
    def visualize(self, obj, **kwargs):
        self._visualize(obj, **kwargs)[0].show()
