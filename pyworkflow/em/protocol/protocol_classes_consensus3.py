# **************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.object import Float, Object, String, Float, Integer
from pyworkflow.utils.path import cleanPath
from pyworkflow.protocol.params import MultiPointerParam

from pyworkflow.em import EMSet, Class3D
from protocol import EMProtocol

# For the viewer part
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView


class ProtClassesConsensus3(EMProtocol):
    """ Compare several SetOfClasses
    """
    _label = 'consensus 3Dclasses (new)'

    interDB = []

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputMultiClasses', MultiPointerParam, important=True,
                      label="Input Classes", pointerClass='SetOfClasses',
                      help='Select several sets of classes where '
                           'to evaluate the consensus.')

        form.addParallelSection(threads=4, mpi=0)

#--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Inserting one step for each couple of setOfClasses comparision
        """

        self.setsLengths = []
        combiList = []

        indx1 = 0
        for pointer1 in self.inputMultiClasses:
            indx1 += 1
            set1 = pointer1.get()

            idSet1 = set1.getObjId()
            
            self.setsLengths.append(set1.getSize())
            indx2 = 0
            for pointer2 in self.inputMultiClasses:
                indx2 += 1
                set2 = pointer2.get()

                idSet2 = set2.getObjId()
                
                idsTup = (indx1,indx2)

                if not (idsTup in combiList or idsTup[::-1] in combiList or 
                        idsTup[0]==idsTup[1]):
                    
                    combiList.append(idsTup)

                    self._insertFunctionStep('compareClassesStep', indx1, indx2)

                
            
            

        self._insertFunctionStep('createOutputStep')

    def compareClassesStep(self, idxSetCls1, idxSetCls2):

        print('idxSetCls1 = %d'%(idxSetCls1))
        print('idxSetCls2 = %d'%(idxSetCls2))
        
        set1 = self.inputMultiClasses[idxSetCls1-1].get()
        set2 = self.inputMultiClasses[idxSetCls2-1].get()

        inx = 0
        for cls1 in set1:
            ids1 = cls1.getIdSet()
            cl1Id = cls1.getObjId()
            ind1Tuple = (idxSetCls1, cl1Id)
            for cls2 in set2:
                ids2 = cls2.getIdSet()

                inter = ids1.intersection(ids2)
                interPop = len(inter)

                if interPop/len(ids1)>interPop/len(ids2):
                    repCls = cls1
                else:
                    repCls = cls2

                cl2Id = cls2.getObjId()
                ind2Tuple = (idxSetCls2, cl2Id)
                interTuple = (len(inter), inter, ind1Tuple, ind2Tuple, repCls)
                self.interDB.append(interTuple)
                inx += 1


    def createOutputStep(self):

        for iii in range(0,len(self.interDB)):

            printTuple = (self.interDB[iii][0], self.interDB[iii][2], 
                          self.interDB[iii][3], self.interDB[iii][4])
            print(printTuple)

        print(len(self.interDB))


        self.interDB.sort(key=lambda e: e[0], reverse=True)

        for iii in range(0,len(self.interDB)):

            printTuple = (self.interDB[iii][0], self.interDB[iii][2], 
                          self.interDB[iii][3], self.interDB[iii][4])
            print(printTuple)



        outputFn = self._getPath('consensus.sqlite')
        cleanPath(outputFn)
        outputSet = EMSet(filename=outputFn)

        for clInx in range(0, min(self.setsLengths)):

            partIds = coinDict[clInx][2]

            cl = Class3D()
            cl.copyInfo(coinDict[clInx][3])
            cl.setRepresentative(coinDict[clInx][3].getRepresentative())

            for i, elem in enumerate(coinDict[clInx][3]):
                if i in partIds:
                    cl.append(elem)


            outputSet.append(cl)

                
        self._defineOutputs(outputConsensus=outputSet)

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors
    
    
class ViewerClassesConsensus(Viewer):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtClassesConsensus3]
    
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
