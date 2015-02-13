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

import os
from pyworkflow.utils.utils import envVarOn

from pyworkflow.manager import Manager
from pyworkflow.config import MenuConfig, ProjectSettings, SettingList
from pyworkflow.project import Project
from pyworkflow.gui import Message
from pyworkflow.gui.browser import FileBrowserWindow
from pyworkflow.gui.plotter import Plotter
import threading

import shlex
from pyworkflow.gui.text import _open_cmd
import SocketServer
import matplotlib.pyplot as plt
# Import possible Object commands to be handled
from pyworkflow.em.showj import OBJCMD_NMA_PLOTDIST, OBJCMD_NMA_VMD, OBJCMD_MOVIE_ALIGNPOLAR, OBJCMD_MOVIE_ALIGNCARTESIAN, OBJCMD_MOVIE_ALIGNPOLARCARTESIAN


from base import ProjectBaseWindow, VIEW_PROTOCOLS, VIEW_PROJECTS



class ProjectWindow(ProjectBaseWindow):
    """ Main window for working in a Project. """
    def __init__(self, path, master=None):
        # Load global configuration
        self.projName = Message.LABEL_PROJECT + os.path.basename(path)
        self.projPath = path
        self.loadProject()

        # TODO: put the menu part more nicely. From here:
        self.menuList = SettingList()
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

        self.menuList.append(menu)
        self.menuCfg = self.menuList.getItem()
        # TODO: up to here

        self.icon = self.generalCfg.icon.get()
        self.selectedProtocol = None
        self.showGraph = False
        Plotter.setBackend('TkAgg')
        ProjectBaseWindow.__init__(self, self.projName, master, icon=self.icon, minsize=(900,500))
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
        self.protCfg = self.settings.getCurrentProtocolMenu()

    #
    # The next functions are callbacks from the menu options.
    # See how it is done in pyworkflow/gui/gui.py:Window._addMenuChilds()
    #
    def onBrowseFiles(self):
        # Project -> Browse files
        FileBrowserWindow("Browse Project files",
                          self, self.project.getPath(''), 
                          selectButton=None  # we are not going to select nothing
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

    def plotSqlite(self, argv):
        dbName = argv[0]
        dbPreffix = argv[1]
        columns = argv[2].split()
        colors = argv[3].split()
        lines = argv[4].split()
        markers = argv[5].split()
        xcolumn = argv[6]
        ylabel = argv[7]
        xlabel = argv[8]
        title = argv[9]
        bins = argv[10]

        from pyworkflow.mapper.sqlite import SqliteFlatDb
        db = SqliteFlatDb(dbName=dbName, tablePrefix=dbPreffix)
        setClassName = db.getProperty('self') # get the set class name
        from pyworkflow.em import getObjects
        setObj = getObjects()[setClassName](filename=dbName, prefix=dbPreffix)

        if xcolumn:
            xvalues = []

            for obj in setObj:
                if hasattr(obj, xcolumn):
                    value = getattr(obj, xcolumn)
                elif xcolumn == 'id':
                    id = int(obj.getObjId())
                    xvalues.append(id)

        else:
            xvalues = range(0, setObj.getSize())


        i = 0
        for column in columns:
            yvalues = []

            for obj in setObj.iterItems(orderBy=column):
                if hasattr(obj, column):

                    value = getattr(obj, column)
                    yvalues.append(value.get())
            color = colors[i]
            line = lines[i]
            if bins:
                #self.queue.put(lambda: getPlotter(xlabel, ylabel, title, xvalues, yvalues, color, line, marker=None, bins=int(bins)).show())
                getPlotter(xlabel, ylabel, title, xvalues, yvalues, color, line, marker=None, bins=int(bins)).show()
            else:
                marker = (markers[i] if not markers[i] == 'none' else None)
                self.queue.put(lambda: getPlotter(xlabel, ylabel, title, xvalues, yvalues, color, line, marker=marker).show())
            i += 1

    def runObjectCommand(self, args):
        from pyworkflow.em.packages.xmipp3.nma.viewer_nma import createDistanceProfilePlot
        from pyworkflow.em.packages.xmipp3.protocol_movie_alignment import createPlots, PLOT_POLAR, PLOT_CART, PLOT_POLARCART
        from pyworkflow.em.packages.xmipp3.nma.viewer_nma import createVmdView
        cmd = args[0]
        protocolId = int(args[1])
        objId = int(args[2])
        project = self.project
        protocol = project.mapper.selectById(protocolId)

        #Plotter.setBackend('TkAgg')
        print cmd
        if cmd == OBJCMD_NMA_PLOTDIST:
            self.queue.put(lambda: createDistanceProfilePlot(protocol, modeNumber=objId).show())

        elif cmd == OBJCMD_NMA_VMD:
            vmd = createVmdView(protocol, modeNumber=objId)
            vmd.show()

        elif cmd == OBJCMD_MOVIE_ALIGNPOLAR:
            self.queue.put(lambda: createPlots(PLOT_POLAR, protocol, objId))

        elif cmd == OBJCMD_MOVIE_ALIGNCARTESIAN:
            self.queue.put(lambda: createPlots(PLOT_CART, protocol, objId))

        elif cmd == OBJCMD_MOVIE_ALIGNPOLARCARTESIAN:
            self.queue.put(lambda: createPlots(PLOT_POLARCART, protocol, objId))


def getPlotter(xlabel, ylabel, title, xvalues, yvalues, color, line, marker=None, bins=None):
    print 'get plotter'
    plotter = Plotter(windowTitle=title)
    a = plotter.createSubPlot(title, xlabel, ylabel)

    if bins:
       a.hist(yvalues, bins=int(bins), color=color, linestyle=line)
    else:
        a.plot(xvalues, yvalues, color, marker=marker, linestyle=line)
    return plotter


class ProjectManagerWindow(ProjectBaseWindow):
    """ Windows to manage all projects. """
    def __init__(self, **args):
        # Load global configuration
        settings = ProjectSettings()
        settings.loadConfig()

        # TODO: put the menu part more nicely. From here:
        menuList = SettingList()
        menu = MenuConfig()

        confMenu = menu.addSubMenu('Configuration')
        confMenu.addSubMenu('General', 'general')
        confMenu.addSubMenu('Hosts', 'hosts')
        confMenu.addSubMenu('Protocols', 'protocols')

        helpMenu = menu.addSubMenu('Help')
        helpMenu.addSubMenu('Online help', 'online_help', icon='fa-external-link.png')
        helpMenu.addSubMenu('About', 'about', icon='fa-question-circle.png')

        menuList.append(menu)
        # TODO: up to here

        self.menuCfg = menuList.getItem()
        self.generalCfg = settings.getConfig()
        
        ProjectBaseWindow.__init__(self, Message.LABEL_PROJECTS, minsize=(750, 500), **args)
        self.manager = Manager()
        
        self.switchView(VIEW_PROJECTS)

    #
    # The next functions are callbacks from the menu options.
    # See how it is done in pyworkflow/gui/gui.py:Window._addMenuChilds()
    #
    def onGeneral(self):
        # Config -> General
        _open_cmd('%s/.config/scipion/scipion.conf' % os.environ['HOME'])

    def onHosts(self):
        # Config -> Hosts
        _open_cmd('%s/.config/scipion/hosts.conf' % os.environ['HOME'])

    def onProtocols(self):
        # Config -> Protocols
        _open_cmd('%s/.config/scipion/protocols.conf' % os.environ['HOME'])



class ProjectTCPServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):
    pass



class ProjectTCPRequestHandler(SocketServer.BaseRequestHandler):

    def handle(self):
        project = self.server.project
        window = self.server.window
        msg = self.request.recv(1024)
        tokens = shlex.split(msg)
        print msg
        if msg.startswith('run protocol'):
            protocolName = tokens[2]
            from pyworkflow.em import getProtocols
            protocolClass = getProtocols()[protocolName]
            # Create the new protocol instance and set the input values
            protocol = project.newProtocol(protocolClass)

            for token in tokens[3:]:
                print token
                param, value = token.split('=')
                attr = getattr(protocol, param, None)
                if param == 'label':
                    protocol.setObjLabel(value)
                elif attr.isPointer():
                    obj = project.getObject(int(value))
                    attr.set(obj)
                elif value:
                    attr.set(value)
            project.launchProtocol(protocol, wait=True)
        if msg.startswith('run function'):
            functionName = tokens[2]
            functionPointer = getattr(window, functionName)
            functionPointer(tokens[3:])
        else:
            answer = 'no answer available'
            self.request.sendall(answer + '\n')

