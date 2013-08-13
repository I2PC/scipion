#!/usr/bin/env python
# To run only the tests in this file, use:
# python -m unittest test_mappers
# To run a single test,
# python -m unittest test_mappers.TestMappers.test_connectUsing

import os
import os.path
import unittest
import pyworkflow.mapper.postgresql
from pyworkflow.object import *


# @see test_object.TestPyworkflow.test_SqliteMapper
class TestMappers(unittest.TestCase):

    def getScipionHome(self):
        if "SCIPION_HOME" not in os.environ:
            raise Exception("SCIPION_HOME is not defined as environment variable")

        return os.environ["SCIPION_HOME"]


    def test_PostgresqlMapper(self):
        # Note: general-purpose exception-handling is handled by Pyunit
        # psql conn parameters? psql python library? !!!
        mapper = pyworkflow.mapper.postgresql.PostgresqlMapper()

        i = Integer(4)

        mapper.insert(i)

        mapper.commit()

        objects = mapper.selectAll()
        self.assertEqual(len(objects),1)

    def test_connectUsing(self):
        db= pyworkflow.mapper.postgresql.PostgresqlDb()
        dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
        db.connectUsing(dbconfig)
