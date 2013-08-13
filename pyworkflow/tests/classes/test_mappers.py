#!/usr/bin/env python
# To run only the tests in this file, use:
# python -m unittest test_mappers

import unittest
from pyworkflow.mapper import *


# @see test_object.TestPyworkflow.test_SqliteMapper
class TestMappers(unittest.TestCase):

    def test_PostgresqlMapper(self):
        # Note: general-purpose exception-handling is handled by Pyunit
        # psql conn parameters? psql python library? !!!
        mapper = PostgresqlMapper()

        i = Integer(4)
        mapper.insert(i)
        mapper.commit()
