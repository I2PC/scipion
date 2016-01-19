#!/usr/bin/env python
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
"""
Main Project window implementation.
It is composed by three panels:
1. Left: protocol tree.
2. Right upper: VIEWS (Data/Protocols)
3. Summary/Details
"""

import os, sys
import threading
import shlex
from pyworkflow.utils import envVarOn, getLocalHostName, getLocalUserName
from pyworkflow.manager import Manager
from pyworkflow.config import MenuConfig, ProjectSettings
from pyworkflow.project import Project
from pyworkflow.gui import Message
from pyworkflow.gui.browser import FileBrowserWindow
from pyworkflow.em.plotter import plotFile
from pyworkflow.gui.plotter import Plotter
from pyworkflow.gui.text import _open_cmd
import SocketServer

# Import possible Object commands to be handled
from pyworkflow.em.showj import (OBJCMD_NMA_PLOTDIST, OBJCMD_NMA_VMD, 
                                 OBJCMD_MOVIE_ALIGNCARTESIAN, OBJCMD_CTFFIND4,
                                 OBJCMD_GCTF)
from base import ProjectBaseWindow, VIEW_PROTOCOLS, VIEW_PROJECTS



class ProjectWindow(ProjectBaseWindow):
    """ Main window for working in a Project. """
    def __init__(self, path, master=None):
        # Load global configuration
        self.projName = Message.LABEL_PROJECT + os.path.basename(path)
        try:
            projTitle = '%s (%s on %s)' % (self.projName, 
                                           getLocalUserName(), 
                                           getLocalHostName())
        except Exception:
            projTitle = self.projName 
            
        self.projPath = path
        self.loadProject()

        # TODO: put the menu part more nicely. From here:
        menu = MenuConfig()

        projMenu = menu.addSubMenu('Project')
        projMenu.addSubMenu('Browse files', 'browse', icon='fa-folder-open.png')
        projMenu.addSubMenu('Remove temporary files', 'delete', icon='fa-trash-o.png')
        projMenu.addSubMenu('', '') # add separator
        projMenu.addSubMenu('Import workflow', 'load_workflow', icon='fa-download.png')
        projMenu.addSubMenu('Export tree graph', 'export_tree')
        projMenu.addSubMenu('', '') # add separator
        projMenu.addSubMenu('Exit', 'exit', icon='fa-sign-out.png')

        helpMenu = menu.addSubMenu('Help')
        helpMenu.addSubMenu('Online help', 'online_help', icon='fa-external-link.png')
        helpMenu.addSubMenu('About', 'about', icon='fa-question-circle.png')

        self.menuCfg = menu
        # TODO: up to here

        self.icon = self.generalCfg.icon.get()
        self.selectedProtocol = None
        self.showGraph = False
        Plotter.setBackend('TkAgg')
        ProjectBaseWindow.__init__(self, projTitle, master, icon=self.icon, minsize=(900,500))
        self.switchView(VIEW_PROTOCOLS)

        self.initProjectTCPServer()#Socket thread to communicate with clients


    def createHeaderFrame(self, parent):
        """Create the header and add the view selection frame at the right."""
        header = ProjectBaseWindow.createHeaderFrame(self, parent)
        self.addViewList(header)
        return header

    def getSettings(self):
        return self.settings
    
    def saveSettings(self):
        self.settings.write()
        
    def _onClosing(self):
        try:
            self.saveSettings()
        except Exception, ex:
            print Message.NO_SAVE_SETTINGS + str(ex) 
        ProjectBaseWindow._onClosing(self)
     
    def loadProject(self):
        self.project = Project(self.projPath)
        self.project.load()
        self.settings = self.project.getSettings()
        self.generalCfg = self.settings.getConfig()
        self.protCfg = self.project.getCurrentProtocolView()

    #
    # The next functions are callbacks from the menu options.
    # See how it is done in pyworkflow/gui/gui.py:Window._addMenuChilds()
    #
    def onBrowseFiles(self):
        # Project -> Browse files
        FileBrowserWindow("Browse Project files",
                          self, self.project.getPath(''), 
                          selectButton=None  # we will select nothing
                          ).show()

    def onRemoveTemporaryFiles(self):
        # Project -> Remove temporary files
        tmpPath = os.path.join(self.project.path, self.project.tmpPath)
        n = 0
        try:
            for fname in os.listdir(tmpPath):
                fpath = "%s/%s" % (tmpPath, fname)
                if os.path.isfile(fpath):
                    os.remove(fpath)
                    n += 1
                # TODO: think what to do with directories. Delete? Report?
            self.showInfo("Deleted content of %s -- %d file(s)." % (tmpPath, n))
        except Exception as e:
            self.showError(str(e))
        
    def _loadWorkflow(self, obj):
        try:
            self.project.loadProtocols(obj.getPath())
            self.getViewWidget().updateRunsGraph(True, reorganize=True)
        except Exception, ex:
            self.showError(str(ex))
            
    def onImportWorkflow(self):
        FileBrowserWindow("Select workflow .json file",
                          self, self.project.getPath(''),
                          onSelect=self._loadWorkflow,
                          selectButton='Import'
                          ).show()

    def onExportTreeGraph(self):
        runsGraph = self.project.getRunsGraph(refresh=True)
        useId = not envVarOn('SCIPION_TREE_NAME')
        runsGraph.printDot(useId=useId)
        if useId:
            print "\nexport SCIPION_TREE_NAME=1 # to use names instead of ids"
        else:
            print "\nexport SCIPION_TREE_NAME=0 # to use ids instead of names"

    def initProjectTCPServer(self):
        server = ProjectTCPServer((self.project.address, self.project.port), ProjectTCPRequestHandler)
        server.project = self.project
        server.window = self
        server_thread = threading.Thread(target=server.serve_forever)
        # Exit the server thread when the main thread terminates
        server_thread.daemon = True
        server_thread.start()

    def schedulePlot(self, path, *args):
        self.enqueue(lambda: plotFile(path, *args).show())    

    def runObjectCommand(self, cmd, inputStrId, objStrId):
        try:
            
            from pyworkflow.em.packages.xmipp3.nma.viewer_nma import createDistanceProfilePlot
            from pyworkflow.em.packages.xmipp3.protocol_movie_alignment import createPlots 
            from pyworkflow.em.protocol.protocol_movies import PLOT_CART
            from pyworkflow.em.packages.xmipp3.nma.viewer_nma import createVmdView
            objId = int(objStrId)
            project = self.project
            if os.path.isfile(inputStrId) and os.path.exists(inputStrId):
                from pyworkflow.em import loadSetFromDb
                inputObj = loadSetFromDb(inputStrId)
            else:
                inputId = int(inputStrId)
                inputObj = project.mapper.selectById(inputId)
            #Plotter.setBackend('TkAgg')
            if cmd == OBJCMD_NMA_PLOTDIST:
                self.enqueue(lambda: createDistanceProfilePlot(inputObj, modeNumber=objId).show())
    
            elif cmd == OBJCMD_NMA_VMD:
                vmd = createVmdView(inputObj, modeNumber=objId)
                vmd.show()
    
    #         elif cmd == OBJCMD_MOVIE_ALIGNPOLAR:
    #             self.enqueue(lambda: createPlots(PLOT_POLAR, inputObj, objId))
    
            elif cmd == OBJCMD_MOVIE_ALIGNCARTESIAN:
                self.enqueue(lambda: createPlots(PLOT_CART, inputObj, objId))
    
            #elif cmd == OBJCMD_MOVIE_ALIGNPOLARCARTESIAN:
            #    self.enqueue(lambda: createPlots(PLOT_POLARCART, inputObj, objId))
            
            elif cmd == OBJCMD_CTFFIND4:
                from pyworkflow.em.packages.grigoriefflab.viewer import createCtfPlot
                self.enqueue(lambda: createCtfPlot(inputObj, objId))
                
            elif cmd == OBJCMD_GCTF:
                from pyworkflow.em.packages.gctf.viewer import createCtfPlot

        except Exception, ex:
            print "There was an error executing object command !!!:"
            print  ex
    



    def recalculateCTF(self, inputObjId, sqliteFile):
        """ Load the project and launch the protocol to
        create the subset.
        """
        # Retrieve project, input protocol and object from db
        project = self.project
        inputObj = project.mapper.selectById(int(inputObjId))
        parentProtId = inputObj.getObjParentId()
        parentProt = project.mapper.selectById(parentProtId)
        protDep = project._getProtocolsDependencies([parentProt])
        if protDep:
            prot = project.copyProtocol(parentProt)
            prot.continueRun.set(parentProt)
        else:
            prot = parentProt
            prot.isFirstTime.set(True)
        
        # Define the input params of the new protocol
        prot.recalculate.set(True)
        prot.sqliteFile.set(sqliteFile)
        # Launch the protocol
        self.getViewWidget().executeProtocol(prot)


class ProjectManagerWindow(ProjectBaseWindow):
    """ Windows to manage all projects. """
    def __init__(self, **args):
        # Load global configuration
        settings = ProjectSettings()

        # TODO: put the menu part more nicely. From here:
        menu = MenuConfig()

        fileMenu = menu.addSubMenu('File')
        fileMenu.addSubMenu('Browse files', 'browse', icon='fa-folder-open.png')
        fileMenu.addSubMenu('Exit', 'exit', icon='fa-sign-out.png')

        confMenu = menu.addSubMenu('Configuration')
        confMenu.addSubMenu('General', 'general')
        confMenu.addSubMenu('Hosts', 'hosts')
        confMenu.addSubMenu('Protocols', 'protocols')

        helpMenu = menu.addSubMenu('Help')
        helpMenu.addSubMenu('Online help', 'online_help', icon='fa-external-link.png')
        helpMenu.addSubMenu('About', 'about', icon='fa-question-circle.png')

        self.menuCfg = menu
        self.generalCfg = settings.getConfig()

        try:
            title = '%s (%s on %s)' % (Message.LABEL_PROJECTS, 
                                       getLocalUserName(), 
                                       getLocalHostName())
        except Exception:
            title = Message.LABEL_PROJECTS
        
        ProjectBaseWindow.__init__(self, title, minsize=(750, 500), **args)
        self.manager = Manager()
        
        self.switchView(VIEW_PROJECTS)

    #
    # The next functions are callbacks from the menu options.
    # See how it is done in pyworkflow/gui/gui.py:Window._addMenuChilds()
    #
    def onBrowseFiles(self):
        # File -> Browse files
        FileBrowserWindow("Browse files", self,
                          os.environ['SCIPION_USER_DATA'], selectButton=None).show()

    def onGeneral(self):
        # Config -> General
        _open_cmd('%s/.config/scipion/scipion.conf' % os.environ['HOME'])
        _open_cmd('%s/config/scipion.conf' % os.environ['SCIPION_HOME'])

    def onHosts(self):
        # Config -> Hosts
        _open_cmd('%s/config//hosts.conf' % os.environ['SCIPION_HOME'])
        if os.path.exists('%s/.config/scipion/hosts.conf' % os.environ['HOME']):
            _open_cmd('%s/.config/scipion/hosts.conf' % os.environ['HOME'])

    def onProtocols(self):
        # Config -> Protocols
        _open_cmd('%s/config/protocols.conf' % os.environ['SCIPION_HOME'])
        if os.path.exists('%s/.config/scipion/protocols.conf' % os.environ['HOME']):
            _open_cmd('%s/.config/scipion/protocols.conf' % os.environ['HOME'])



class ProjectTCPServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):
    pass



class ProjectTCPRequestHandler(SocketServer.BaseRequestHandler):

    def handle(self):
        try:
            project = self.server.project
            window = self.server.window
            msg = self.request.recv(1024)
            tokens = shlex.split(msg)
            if msg.startswith('run protocol'):
                protocolName = tokens[2]
                from pyworkflow.em import getProtocols
                protocolClass = getProtocols()[protocolName]
                # Create the new protocol instance and set the input values
                protocol = project.newProtocol(protocolClass)

                for token in tokens[3:]:
                    #print token
                    param, value = token.split('=')
                    attr = getattr(protocol, param, None)
                    if param == 'label':
                        protocol.setObjLabel(value)
                    elif attr.isPointer():
                        obj = project.getObject(int(value))
                        attr.set(obj)
                    elif value:
                        attr.set(value)
                #project.launchProtocol(protocol)
                # We need to enqueue the action of execute a new protocol
                # to be run in the same GUI thread and avoid concurrent
                # access to the project sqlite database
                window.getViewWidget().executeProtocol(protocol)
            if msg.startswith('run function'):
                functionName = tokens[2]
                functionPointer = getattr(window, functionName)
                functionPointer(*tokens[3:])
            else:
                answer = 'no answer available'
                self.request.sendall(answer + '\n')
        except:
            import traceback
            traceback.print_stack()

