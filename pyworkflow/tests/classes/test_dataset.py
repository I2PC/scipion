#!/usr/bin/env python

import os, time
import unittest
import pyworkflow.dataset as ds
from pyworkflow.tests import getInputPath



class TestDataSet(unittest.TestCase):
    """ Some tests for DataSet implementation. """

    def setUp(self):
        pass
        
    def test_Table(self):
        table = ds.Table(ds.Column('x', int, 5),
                         ds.Column('y', float, 0.0),
                         ds.Column('name', str))
        
        # Add a row to the table
        table.addRow(1, x=12, y=11.0, name='jose')
        table.addRow(2, x=22, y=21.0, name='juan')
        table.addRow(3, x=32, y=31.0, name='pedro')
        # Expect an exception, since name is not provided and have not default
        self.assertRaises(Exception, table.addRow, 100, y=3.0)
        
        row = table.getRow(1)
        print row
        
        self.assertEqual(table.getSize(), 3, "Bad table size")
        
        # Update a value of a row
        table.updateRow(1, name='pepe')        
        row = table.getRow(1)
        print row
        self.assertEqual(row.name, 'pepe', "Error updating name in row")

        print "Table:"
        print table

#    def test_XmippDataSet(self):
#        """ Create a table from a metadata. """
#        from pyworkflow.em.packages.xmipp3 import XmippDataSet
#        import xmipp
#        mdPath = getInputPath('showj', 'tux_vol.xmd')
#
#        xds = XmippDataSet(mdPath)
#        
#        tableNames = xds.listTables()
#        print '\ntables: ', tableNames
#        
#        tn = tableNames[0]
#        
#        table = xds.getTable(tn)
#        
#        print "\nTable '%s':" % tn
#        print table
#        
#        md = xds._convertTableToMd(table)
#        print md
        
if __name__ == '__main__':
    unittest.main()
    