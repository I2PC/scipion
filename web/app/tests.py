"""
This file demonstrates writing tests using the unittest module. These will pass
when you run "manage.py test".
"""
import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *

from django.test.client import Client
from django.test.client import RequestFactory
from django.test import TestCase

from pyworkflow.web.pages import settings
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from pyworkflow.utils.utils import prettyDate

from pyworkflow.web.app import views

"""
django.test.TestCase is a subclass of unittest.TestCase that runs each test 
inside a transaction to provide isolation
"""

#
# class SimpleTest(TestCase):
#    def test_basic_addition(self):
#        """
#        Tests that 1 + 1 always equals 2.
#        """
#        self.assertEqual(1 + 1, 2)

class TestWebWorkflow(TestCase):
    
    def setUp(self):
        # Every test needs a client.
        self.client = Client()
           
    def testProjects(self):
        
        # Using the view.py method
        
        req = RequestFactory().get('/projects')
        views.projects(req)
        
        # Using the client instance
        
        manager = Manager()
                        
        projects = manager.listProjects()
        for p in projects:
            p.pTime = prettyDate(p.mTime)
                
        css_path = os.path.join(settings.STATIC_URL, 'css/projects_style.css')
        
        context = {'projects': projects,
                   'css':css_path}
   
        response = self.client.get('/projects/', context)
        
        print response

