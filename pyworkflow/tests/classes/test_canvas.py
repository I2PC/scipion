#!/usr/bin/env python
# To run only the tests in this file, use:
# python -m unittest test_canvas -v


import unittest
from pyworkflow.gui.canvas import Canvas
import Tkinter
import math
from sets import Set

class TestCanvas(unittest.TestCase):

    def distance(self,c1,c2):
        return math.hypot(c2[0] - c1[0], c2[1] - c1[1])

    def allEqual(self,list):
        return not list or list.count(list[0]) == len(list)

    def allDifferent(self,list):
        return not list or len(list) == len(Set(list))

    def test_connectorsCoords(self):
        root = Tkinter.Tk()
        canvas = Canvas(root, width=800, height=600)
        tb1 = canvas.createTextbox("Project", 100, 100, "blue")
        print str(tb1.getDimensions())

        connectorsCoords=canvas.getConnectorsCoordinates(tb1)
        self.assertTrue(self.allDifferent(connectorsCoords))

        print connectorsCoords

        distances={}
        for i in range(len(connectorsCoords)-1):
            distances[i]=self.distance(connectorsCoords[i],connectorsCoords[i+1])

        self.assertTrue(self.allEqual(distances.values()))
        self.assertNotEqual(distances[0],0)

