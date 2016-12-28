# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.tests import *
import pyworkflow.em.metadata as md



class TestMetaData(unittest.TestCase):
    
    _labels = [WEEKLY]

    def _newMd(self):
        md0 = md.MetaData()
        n = 5
        xcoor = range(n)
        ycoor = [x*x for x in xcoor]
        for i in range(n):
            self._addRow(md0, '%02d@proj.stk' % i, xcoor[i], ycoor[i])
        return md0

    def _addRow(self, md0, imageFn, xcoor, ycoor):
        objId = md0.addObject()
        md0.setValue(md.MDL_IMAGE, imageFn, objId)
        md0.setValue(md.MDL_XCOOR, xcoor, objId)
        md0.setValue(md.MDL_YCOOR, ycoor, objId)

    def test_removeDuplicates(self):
        md0 = self._newMd()
        md1 = self._newMd()

        # If removing without labels, this metadata should remain the same
        md1.removeDuplicates()
        self.assertEqual(md0, md1)

        # We can use labels for removeDuplicates
        self._addRow(md1, '00@proj.stk', 0, 0)
        md1.removeDuplicates()
        self.assertEqual(md0, md1)

        self._addRow(md1, '06@proj.stk', 0, 0)
        self.assertNotEqual(md0, md1)
        md1.removeDuplicates()
        self.assertNotEqual(md0, md1)
        md1.removeDuplicates(md.MDL_XCOOR)
        self.assertEqual(md0, md1)

        md1.clear()
        self._addRow(md1, '00@proj.stk', 0, 1)
        self._addRow(md1, '00@proj.stk', 0, 2)
        self._addRow(md1, '00@proj.stk', 0, 3)
        md1.removeDuplicates(md.MDL_IMAGE)
        self.assertEqual(md1.size(), 1)
        objId = md1.firstObject()
        self.assertEqual(md1.getValue(md.MDL_YCOOR, objId), 1)

