# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

from __future__ import print_function
import os
import sys
import shlex
import ast
from threading import Thread
from multiprocessing.connection import Client
from numpy import flipud
import socket
try:  # python 2
    import Tkinter as tk
    import tkFont
    import ttk
except ImportError:  # Python 3
    import tkinter as tk
    import tkinter.font as tkFont
    import tkinter.ttk as ttk

import pyworkflow as pw
from pyworkflow.viewer import View, Viewer, CommandView, DESKTOP_TKINTER
from pyworkflow.utils import Environ, runJob
from pyworkflow.utils import getFreePort
from pyworkflow.gui.matplotlib_image import ImageWindow

# From pyworkflow.em level
import showj
import metadata as md
from data import PdbFile
from convert import ImageHandler

import xmipp

# from viewer_fsc import FscViewer
# from viewer_pdf import PDFReportViewer
# from viewer_monitor_summary import ViewerMonitorSummary
# from protocol.monitors.protocol_monitor_ctf import ProtMonitorCTFViewer
# from protocol.monitors.protocol_monitor_system import ProtMonitorSystemViewer

# ------------------------ Some common Views ------------------


class DataView(View):
    """ Wrapper the arguments to showj (either web or desktop). """
    def __init__(self, path, viewParams={}, **kwargs):
        View.__init__(self)
        self._memory = '2g'
        self._loadPath(path)
        self._env = kwargs.get('env', {})
        self._viewParams = viewParams

    def _loadPath(self, path):
        self._tableName = None

        # If path is a tuple, we will convert to the filename format
        # as expected by Showj
        if isinstance(path, tuple):
            self._path = ImageHandler.locationToXmipp(path)
        # Check if there is a table name with @ in path
        # in that case split table name and path
        # table names can never starts with a number
        # this is considering an image inside an stack
        elif isinstance(path, basestring):
            if '@' in path and path[0] not in '0123456789':
                self._tableName, self._path = path.split('@')
            else:
                self._path = path
        else:
            raise Exception("Invalid input path, "
                            "should be 'string' or 'tuple'")

    def show(self):
        showj.runJavaIJapp(self._memory, 'xmipp.viewer.scipion.ScipionViewer',
                           self.getShowJParams(), env=self._env)

    def getShowJParams(self):
        tableName = '%s@' % self._tableName if self._tableName else ''
        params = '-i "%s%s"' % (tableName, self._path)
        for key, value in self._viewParams.items():
            params = "%s --%s %s" % (params, key, value)

        return params

    def getShowJWebParams(self):

        """
        ===========OLD SHOWJ WEB DOCUMENTATION=======================
        Extra parameters can be used to configure table layout and set render
        function for a column
        Default layout configuration is set in ColumnLayoutProperties method in
        layout_configuration.py

        Parameters are formed by: [label]___[property]: [value].
        E.g.: id___visible:True or micrograph___renderFunc:"get_image_psd"
        Properties to be configured are:
            visible: Defines if this column is displayed
            allowSetVisible: Defines if user can change visible property
            (show/hide this column).
            editable: Defines if this column is editable, ie user can change
            field value.
            allowSetEditable: Defines if user can change editable property
            (allow editing this column).
            renderable: Defines if this column is renderizable, ie it renders
            data column using renderFunc
            allowSetRenderable: Defines if user can change renderable property.
            renderFunc: Function to be used when this field is rendered.
            (it has to be inserted in render_column method)
            extraRenderFunc: Any extra parameters needed for rendering.
            Parameters are passed like in a url ie downsample = 2 &
            lowPass = 3.5

        Example:
        extraParameters["id___visible"]=True
        extraParameters["micrograph___renderFunc"]="get_image_psd"
        extraParameters["micrograph___extraRenderFunc"]="downsample=2"
        ========================================================================
        """
        parameters = {
            showj.MODE,  # FOR MODE TABLE OR GALLERY
            showj.VISIBLE,
            showj.ZOOM,
            showj.ORDER,
            showj.RENDER,
            showj.SORT_BY}

        params = {}

        for key, value in self._viewParams.items():
            if key in parameters:
                if key == 'mode' and value == 'metadata':
                    value = 'table'
                params[key] = value

        return params

    def getPath(self):
        return self._path

    def getTableName(self):
        return self._tableName


class ObjectView(DataView):
    """ Wrapper to DataView but for displaying Scipion objects. """
    def __init__(self, project, inputid, path, other='', viewParams={},
                 **kwargs):
        DataView.__init__(self, path, viewParams, **kwargs)
        self.type = type
        self.port = project.port
        self.inputid = inputid
        self.other = other

    def getShowJParams(self):
        # Add the scipion parameters over the normal showj params
        return '%s --scipion %s %s %s' % (DataView.getShowJParams(self),
                                          self.port, self.inputid, self.other)

    def show(self):
        showj.runJavaIJapp(self._memory, 'xmipp.viewer.scipion.ScipionViewer',
                           self.getShowJParams(), env=self._env)


class MicrographsView(ObjectView):
    """ Customized ObjectView for SetOfCTF objects . """
    # All extra labels that we want to show if present in the CTF results
    RENDER_LABELS = ['thumbnail._filename', 'psdCorr._filename',
                     'plotGlobal._filename']
    EXTRA_LABELS = ['_filename']

    def __init__(self, project, micSet, other='', **kwargs):
        first = micSet.getFirstItem()

        first.printAll()

        def existingLabels(labelList):
            print("labelList: ", labelList)
            return ' '.join([l for l in labelList if first.hasAttributeExt(l)])

        renderLabels = existingLabels(self.RENDER_LABELS)
        extraLabels = existingLabels(self.EXTRA_LABELS)
        labels = 'id enabled %s %s' % (renderLabels, extraLabels)

        viewParams = {showj.MODE: showj.MODE_MD,
                      showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.ZOOM: 50
                      }

        if renderLabels:
            viewParams[showj.RENDER] = renderLabels

        inputId = micSet.getObjId() or micSet.getFileName()
        ObjectView.__init__(self, project,
                            inputId, micSet.getFileName(), other,
                            viewParams, **kwargs)


class CtfView(ObjectView):
    """ Customized ObjectView for SetOfCTF objects . """
    # All extra labels that we want to show if present in the CTF results
    PSD_LABELS = ['_micObj.thumbnail._filename', '_psdFile',
                  '_xmipp_enhanced_psd', '_xmipp_ctfmodel_quadrant',
                  '_xmipp_ctfmodel_halfplane', '_micObj.plotGlobal._filename']
    EXTRA_LABELS = ['_ctffind4_ctfResolution', '_xmipp_ctfCritFirstZero',
                    ' _xmipp_ctfCritCorr13', '_xmipp_ctfCritFitting',
                    '_xmipp_ctfCritNonAstigmaticValidity',
                    '_xmipp_ctfCritCtfMargin', '_xmipp_ctfCritMaxFreq']

    def __init__(self, project, ctfSet, other='', **kwargs):
        first = ctfSet.getFirstItem()

        def existingLabels(labelList):
            return ' '.join([l for l in labelList if first.hasAttributeExt(l)])

        psdLabels = existingLabels(self.PSD_LABELS)
        extraLabels = existingLabels(self.EXTRA_LABELS)
        labels = 'id enabled %s _defocusU _defocusV ' % psdLabels
        labels += '_defocusAngle _defocusRatio %s  _micObj._filename' % \
                  extraLabels

        viewParams = {showj.MODE: showj.MODE_MD,
                      showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.ZOOM: 50}

        if psdLabels:
            viewParams[showj.RENDER] = psdLabels

        if ctfSet.isStreamOpen():
            viewParams['dont_recalc_ctf'] = ''

        if first.hasAttribute('_ctffind4_ctfResolution'):
            import pyworkflow.em.packages.grigoriefflab.viewer as gviewer
            viewParams[showj.OBJCMDS] = "'%s'" % gviewer.OBJCMD_CTFFIND4

        inputId = ctfSet.getObjId() or ctfSet.getFileName()
        ObjectView.__init__(self, project,
                            inputId, ctfSet.getFileName(), other,
                            viewParams, **kwargs)


class ClassesView(ObjectView):
    """ Customized ObjectView for SetOfClasses. """
    def __init__(self, project, inputid, path, other='',
                 viewParams={}, **kwargs):
        labels = 'enabled id _size _representative._filename'
        defaultViewParams = {showj.ORDER: labels,
                             showj.VISIBLE: labels,
                             showj.RENDER: '_representative._filename',
                             showj.SORT_BY: '_size desc',
                             showj.LABELS: 'id _size',
                             }
        defaultViewParams.update(viewParams)
        ObjectView.__init__(self, project, inputid, path, other,
                            defaultViewParams, **kwargs)


class Classes3DView(ClassesView):
    """ Customized ObjectView for SetOfClasses. """
    def __init__(self, project, inputid, path, other='',
                 viewParams={}, **kwargs):
        defaultViewParams = {showj.ZOOM: '99',
                             showj.MODE: 'metadata'}
        defaultViewParams.update(viewParams)
        ClassesView.__init__(self, project, inputid, path, other,
                             defaultViewParams, **kwargs)


class CoordinatesObjectView(DataView):
    """ Wrapper to View but for displaying Scipion objects. """
    def __init__(self, project, path, outputdir, protocol, pickerProps=None,
                 inTmpFolder=False, **kwargs):
        DataView.__init__(self, path, **kwargs)
        self.project = project
        self.outputdir = outputdir
        self.protocol = protocol
        self.pickerProps = pickerProps
        self.inTmpFolder = inTmpFolder

    def show(self):
        return showj.launchSupervisedPickerGUI(self._path, self.outputdir,
                                               self.protocol,
                                               pickerProps=self.pickerProps,
                                               inTmpFolder=self.inTmpFolder)


class ImageView(View):
    """ Customized ObjectView for SetOfClasses. """
    def __init__(self, imagePath, **kwargs):
        View.__init__(self)
        self._imagePath = os.path.abspath(imagePath)

    def getImagePath(self):
        return self._imagePath
# TODO: delete class TextFileView
        '''
class TextFileView(View):

    def __init__(self, path, tkRoot):
        self.path = path
        self.tkRoot=tkRoot#message box will be painted ABOVE this window

    def show(self):
        """Show text file in default editor, If file does not exists return
        error message"""
        if not os.path.isfile(self.path):
            tkMessageBox.showerror("Refamc Viewer Error",#bar title
                                   "refmac log file not found\n(%s)"
                                   % self.path,#message
                                   parent=self.tkRoot)
            return
        editor = os.getenv('EDITOR')
        if editor:
            os.system(editor + ' ' + self.path)
        else:
            webbrowser.open(self.path)
'''


class TableView(View):
    """ show table, pass values as:
        headerList = ['name', 'surname']
        dataList = [
        ('John', 'Smith') ,
        ('Larry', 'Black') ,
        ('Walter', 'White') ,
        ('Fred', 'Becker')
        ].
        msg = message to be shown at the table top
        title= window title
        height: Specifies the number of rows which should be visible
        width: minimum width in pixels
        fontSize= font size
        padding: cell extra width
        ---------------------
        Alternative way to create a table using showj
        views = []
        labels = '_1 _2'
        emSet = EMSet(filename="/tmp/kk.sqlite")
        emObject = EMObject()
        emObject._1 = String('first parameter')
        emObject._2 = Float(12.)
        emSet.append(emObject)
        emObject = EMObject()
        emObject._1 = String('second parameter')
        emObject._2 = Float(22.)
        emSet.append(emObject)
        emSet.write()
        views.append(ObjectView(self._project,
                                self.protocol.strId(),
                                "/tmp/kk.sqlite",
                                viewParams = {MODE: MODE_MD, ORDER: labels,
                                VISIBLE: labels}))
        return views
"""

    def __init__(self, headerList, dataList,
                 mesg=None, title=None,
                 height=10, width=400,
                 fontSize=16, padding=10):
        # get new widget that has as parent the top level window and set title
        win = tk.Toplevel()
        if title:
            win.wm_title(title)

        # frame to place all other widgets
        frame = tk.Frame(win)

        # make font a little bigger
        # TODO: font size should be general
        font = tkFont.Font(family='fixed', size=fontSize)
        font.metrics()
        fontheight = font.metrics()['linespace']
        style = ttk.Style()
        style.configure('Calendar.Treeview', font=font, rowheight=fontheight)

        # create treeview to store multi list with data
        tree = ttk.Treeview(columns=headerList,
                            show="headings", master=win,
                            style='Calendar.Treeview', height=height
                            )

        # define scrollbars to be added
        if len(dataList) > height:
            ysb = ttk.Scrollbar(orient=tk.VERTICAL, command=tree.yview,
                                master=win)
            # xsb = ttk.Scrollbar(orient=tk.HORIZONTAL,
            # command= tree.xview, master=win)
            # add them to three view
            tree.configure(yscroll=ysb.set)  # , xscroll=xsb.set)

        # create rows and columns
        counterRow = 1
        colWidths = []  # list with maximum width per column
        # create headers
        for col in headerList:
            tree.heading(col, text=col.title())
            # save neede width for this cell
            (w, h) = (font.measure(col.title()), font.metrics("linespace"))
            colWidths.append(w)

        # insert other rows
        # tag rows as odd or even so they may have different background colors
        for item in dataList:
            if counterRow % 2:
                tree.insert('', 'end', values=item, tags=('evenrow',))
            else:
                tree.insert('', 'end', values=item, tags=('oddrow',))
            counterRow += 1
            counterCol = 0
            for i in item:
                (w, h) = (font.measure(i), font.metrics("linespace"))
                if colWidths[counterCol] < w:
                    colWidths[counterCol] = w
                counterCol += 1

        # if width less than sum of column widths expand them
        sumColWid = sum(colWidths) + 20
        if sumColWid < width:
            sumColWid = width
            factor = int(width/sumColWid)+1
            colWidths = [i * factor for i in colWidths]

        for col, colWidth in zip(headerList, colWidths):
            tree.column(col, width=colWidth+padding)

        # color by rows
        tree.tag_configure('evenrow', background='white')
        tree.tag_configure('oddrow', background='light grey')

        # message placed at the window top
        msg = ttk.Label(wraplength=sumColWid, justify="left", anchor="n",
                        padding=(10, 2, 10, 6), text=(mesg), master=win,
                        font=font)

        # set mg in grid 0,0
        msg.grid(row=0, column=0)
        # set tree in grid 1,0
        tree.grid(row=1, column=0)
        # set ysg in grid 1 1
        # but only if number of elements is larger than height
        if len(dataList) > height:
            ysb.grid(row=1, column=1, sticky='ns')
            # xsb.grid(row=2, column=0, sticky='ew')

# ------------------------ Some views and  viewers ------------------------


def getChimeraEnviron():
    """ Return the proper environ to launch chimera.
    CHIMERA_HOME variable is read from the ~/.config/scipion.conf file.
    """
    environ = Environ(os.environ)
    environ.set('PATH', os.path.join(os.environ['CHIMERA_HOME'], 'bin'),
                position=Environ.BEGIN)

    if "REMOTE_MESA_LIB" in os.environ:
        environ.set('LD_LIBRARY_PATH', os.environ['REMOTE_MESA_LIB'],
                    position=Environ.BEGIN)
    return environ


class ChimeraView(CommandView):
    """ View for calling an external command. """
    def __init__(self, inputFile, **kwargs):
        CommandView.__init__(self, 'chimera "%s" &' % inputFile,
                             env=getChimeraEnviron(), **kwargs)


class ChimeraClientView(CommandView):
    """ View for calling an external command. """
    def __init__(self, inputFile, **kwargs):
        self._inputFile = inputFile
        self._kwargs = kwargs

    def show(self):
        if self._kwargs.get('showProjection', False):
            ChimeraProjectionClient(self._inputFile, **self._kwargs)
        else:
            ChimeraAngDistClient(self._inputFile, **self._kwargs)


class ChimeraDataView(ChimeraClientView):

    def __init__(self, dataview, vol, viewParams={}, **kwargs):
        self.dataview = dataview
        self.showjPort = getFreePort()
        self.dataview._viewParams[showj.CHIMERA_PORT] = self.showjPort
        self.dataview._viewParams[showj.MODE] = showj.MODE_MD
        self.dataview._viewParams[showj.INVERTY] = ''
        ChimeraClientView.__init__(self, vol.getFileName(),
                                   showProjection=True,
                                   showjPort=self.showjPort,
                                   voxelSize=vol.getSamplingRate())

    def show(self):
        self.dataview.show()
        ChimeraClientView.show(self)


class ChimeraViewer(Viewer):
    """ Wrapper to visualize PDB object with Chimera. """
    _environments = [DESKTOP_TKINTER]
    _targets = [PdbFile]

    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    def visualize(self, obj, **kwargs):
        cls = type(obj)
        if issubclass(cls, PdbFile):
            fn = obj.getFileName()
            if obj.getPseudoAtoms():
                if hasattr(obj, '_chimeraScript'):
                    fn = obj._chimeraScript.get()
            elif obj.hasVolume():
                pass
                #print "YES"
                #print "obj.hasVolume(): ", obj.hasVolume()
            ChimeraView(fn).show()
            # FIXME: there is an asymmetry between ProtocolViewer and Viewer
            # for the first, the visualize method return a list of View's
            # (that are shown)
            # for the second, the visualize directly shows the objects.
            # the first approach is better
        else:
            raise Exception('ChimeraViewer.visualize: '
                            'can not visualize class: %s' % obj.getClassName())


class ChimeraClient:

    def openVolumeOnServer(self, volume, sendEnd=True):
        self.send('open_volume', volume)
        if self.voxelSize is not None:
            self.send('voxel_size', self.voxelSize)
        if sendEnd:
            self.client.send('end')

    def __init__(self, volfile, sendEnd=True, **kwargs):
        if volfile.endswith('.mrc'):
            volfile += ':mrc'

        self.kwargs = kwargs
        if volfile is None:
            raise ValueError(volfile)
        if '@' in volfile:
            [index, file] = volfile.split('@')
        else:
            file = volfile
        if ':' in file:
            file = file[0: file.rfind(':')]
        if not os.path.exists(file):
            raise Exception("File %s does not exists" % file)

        self.volfile = volfile
        self.voxelSize = self.kwargs.get('voxelSize', None)
        # ChimeraServer is the basic server. There are other
        # than inherit from it.
        serverName = self.kwargs.get('ChimeraServer', 'ChimeraServer')
        self.address = ''
        self.port = getFreePort()

        serverfile = pw.join('em', 'chimera_server.py')
        command = CommandView("chimera --script '%s %s %s' &" %
                              (serverfile, self.port, serverName),
                              env=getChimeraEnviron(),).show()
        self.authkey = 'test'
        self.client = Client((self.address, self.port), authkey=self.authkey)
        self.initVolumeData()
        # self.openVolumeOnServer(self.vol,sendEnd)
        self.openVolumeOnServer(self.vol)
        self.initListenThread()

    def send(self, cmd, data):
        self.client.send(cmd)
        self.client.send(data)

    def initListenThread(self):
            self.listen_thread = Thread(target=self.listen)
            # self.listen_thread.daemon = True
            self.listen_thread.start()

    def listen(self):

        self.listen = True
        try:
            while self.listen:
                msg = self.client.recv()
                self.answer(msg)
        except EOFError:
            print ('Lost connection to server')
        finally:
            self.exit()

    def exit(self):
        self.client.close()

    def initVolumeData(self):
        self.image = xmipp.Image(self.volfile)
        self.image.convert2DataType(md.DT_DOUBLE)
        self.xdim, self.ydim, self.zdim, self.n = self.image.getDimensions()
        self.vol = self.image.getData()

    def answer(self, msg):
        if msg == 'exit_server':
            self.listen = False


class ChimeraAngDistClient(ChimeraClient):

    def __init__(self, volfile, **kwargs):
        self.angularDistFile = kwargs.get('angularDistFile', None)

        if self.angularDistFile:
            fileName = self.angularDistFile
            if '@' in self.angularDistFile:
                fileName = self.angularDistFile.split("@")[1]
            if not (os.path.exists(fileName)):  # check blockname:
                raise Exception("Path %s does not exists" %
                                self.angularDistFile)
        self.spheresColor = kwargs.get('spheresColor', 'red')
        spheresDistance = kwargs.get('spheresDistance', None)
        spheresMaxRadius = kwargs.get('spheresMaxRadius', None)
        ChimeraClient.__init__(self, volfile, sendEnd=False, **kwargs)
        self.spheresDistance = float(spheresDistance) \
            if spheresDistance else 0.75 * max(self.xdim, self.ydim, self.zdim)
        self.spheresMaxRadius = float(spheresMaxRadius) \
            if spheresMaxRadius else 0.02 * self.spheresDistance
        self.loadAngularDist(True)

    def loadAngularDist(self,  sendEnd=True):
        if self.angularDistFile:
            self.readAngularDistFile()
            self.send('command_list', self.angulardist)
        if sendEnd:
            self.client.send('end')

    def readAngularDistFile(self):
        angleRotLabel = md.MDL_ANGLE_ROT
        angleTiltLabel = md.MDL_ANGLE_TILT
        anglePsiLabel = md.MDL_ANGLE_PSI
        mdAngDist = md.MetaData(self.angularDistFile)
        if not mdAngDist.containsLabel(md.MDL_ANGLE_PSI):
            anglePsiLabel = None
            if mdAngDist.containsLabel(md.RLN_ORIENT_PSI):
                angleRotLabel = md.RLN_ORIENT_ROT
                angleTiltLabel = md.RLN_ORIENT_TILT
                anglePsiLabel = md.RLN_ORIENT_PSI

        if not mdAngDist.containsLabel(md.MDL_WEIGHT):
            mdAngDist.fillConstant(md.MDL_WEIGHT, 1.)

        maxweight = mdAngDist.aggregateSingle(md.AGGR_MAX, md.MDL_WEIGHT)
        minweight = mdAngDist.aggregateSingle(md.AGGR_MIN, md.MDL_WEIGHT)
        interval = maxweight - minweight

        self.angulardist = []
        x2 = self.xdim/2
        y2 = self.ydim/2
        z2 = self.zdim/2
        # cofr does not seem to work!
        # self.angulardist.append('cofr %d,%d,%d'%(x2,y2,z2))
        for id in mdAngDist:
            rot = mdAngDist.getValue(angleRotLabel, id)
            tilt = mdAngDist.getValue(angleTiltLabel, id)
            psi = mdAngDist.getValue(anglePsiLabel, id) if anglePsiLabel else 0
            weight = mdAngDist.getValue(md.MDL_WEIGHT, id)
            # Avoid zero division
            weight = 0 if interval == 0 else (weight - minweight)/interval
            weight = weight + 0.5  # add 0.5 to avoid cero weight
            x, y, z = xmipp.Euler_direction(rot, tilt, psi)
            radius = weight * self.spheresMaxRadius

            x = x * self.spheresDistance + x2
            y = y * self.spheresDistance + y2
            z = z * self.spheresDistance + z2
            command = 'shape sphere radius %s center %s,%s,%s color %s ' % \
                      (radius, x, y, z, self.spheresColor)
            self.angulardist.append(command)
            # printCmd(command)


class ChimeraVirusClient(ChimeraClient):

    def __init__(self, volfile, **kwargs):
        self.h = kwargs.get('h', 5)
        if self.h is None:
                self.h = 5
        self.k = kwargs.get('k', 0)
        if self.k is None:
                self.k = 0
        self.sym = kwargs.get('sym', 'i222r')
        if self.sym is None:
                self.sym = 'i222r'
        self.radius = kwargs.get('radius', 100.)
        if self.radius is None:
                self.radius = 100.
        self.spheRadius = kwargs.get('spheRadius', 1.5)
        if self.spheRadius is None:
                self.spheRadius = 1.5
        self.color = kwargs.get('color', 'red')
        if self.color is None:
                self.color = 'red'
        self.linewidth = kwargs.get('linewidth', 1)
        if self.linewidth is None:
                self.linewidth = 1
        self.sphere = kwargs.get('sphere', 0)
        if self.sphere is None:
                self.sphere = 0
        self.shellRadius = kwargs.get('shellRadius', self.spheRadius)
        if self.shellRadius is None:
                self.shellRadius = self.spheRadius
        kwargs['ChimeraServer'] = 'ChimeraVirusServer'
        ChimeraClient.__init__(self, volfile, **kwargs)

    def hkcageCommand(self):
        command = 'hkcage '
        command += '%d %d ' % (self.h, self.k)
        command += 'radius  %f ' % self.radius
        command += 'orientation %s ' % self.sym
        command += 'replace true '
        command += 'color %s ' % self.color
        command += 'linewidth %d ' % self.linewidth
        command += 'sphere %f ' % self.sphere
        return command

    def openVolumeOnServer(self, volume):
        ChimeraClient.openVolumeOnServer(self, volume, sendEnd=False)
        commandList = [self.hkcageCommand()]
        self.send('hk_icosahedron_lattice', (self.h,
                                             self.k,
                                             self.radius,
                                             self.shellRadius,
                                             self.spheRadius,
                                             self.sym,
                                             self.sphere,
                                             self.color))
        # get here the list of vertexes, info will be pass by id later
        self.send('command_list', commandList)
        # endhere
        self.client.send('end')
        return

        # get va with sphere centers

        # va coordinates of  vertex of the triangles inside de canonical
        # triangle
        msg1 = self.client.recv()
        self.va = self.client.recv()
        # self.listToBild(self.va,1.6,msg1+'.bild')
        self.client.send('end')

    # this is an auxiliary file that should be deleted
    def listToBild(self, points, radius, file, color='0 0 1'):
        print ("--------------------------------")
        f = open(file, 'w')
        for point in points:
            print ("\n.color", color, "\n.sphere", point[0]+128.,
                   point[1]+128., point[2]+128., radius, file=f)

    def answer(self, msg):
        ChimeraClient.answer(self, msg)
        if msg == 'motion_stop':
            data = self.client.recv()  # wait for data
            printCmd('reading motion')
            self.motion = data
            printCmd('getting euler angles')
            rot, tilt, psi = xmipp.Euler_matrix2angles(self.motion)
            printCmd('calling rotate')
            self.rotate(rot, tilt, psi)
        elif msg == 'id':
            id = self.client.recv()
            self.listToBild([self.va[id-1]], 2.5, 'id.bild')
        else:
            pass

    def listenShowJ(self):
        """ This function with a very confusing name send messages to
        the chimera server"""

        while True:
            try:
                (clientsocket, address) = self.serversocket.accept()
                # Should be a single message, so no loop
                msg = clientsocket.recv(1024)
                tokens = shlex.split(msg)
                cmd = tokens[0]
                if cmd == 'rotate':
                    rot = float(tokens[1])
                    tilt = float(tokens[2])
                    psi = float(tokens[3])

                    matrix = xmipp.Euler_angles2matrix(rot, tilt, psi)
                elif cmd == 'rotate_matrix':
                    matrixString = tokens[1]
                    matrix = ast.literal_eval(matrixString)
                    matrix = [matrix[0][:3], matrix[1][:3], matrix[2][:3]]
                self.client.send('rotate')
                self.client.send(matrix)
                clientsocket.close()

            except EOFError:
                print ('Lost connection to client')


class ChimeraProjectionClient(ChimeraAngDistClient):

    def __init__(self, volfile, **kwargs):
        print("on chimera projection client")
        ChimeraAngDistClient.__init__(self, volfile, **kwargs)
        self.projection = xmipp.Image()
        self.projection.setDataType(md.DT_DOUBLE)
        # 0.5 ->  Niquiest frequency
        # 2 -> bspline interpolation
        size = self.kwargs.get('size', None)
        defaultSize = self.xdim if self.xdim > 128 else 128
        self.size = size if size else defaultSize
        paddingFactor = self.kwargs.get('paddingFactor', 1)
        maxFreq = self.kwargs.get('maxFreq', 0.5)
        splineDegree = self.kwargs.get('splineDegree', 3)
        self.fourierprojector = \
            xmipp.FourierProjector(self.image, paddingFactor,
                                   maxFreq, splineDegree)
        self.fourierprojector.projectVolume(self.projection, 0, 0, 0)
        self.showjPort = self.kwargs.get('showjPort', None)
        self.iw = ImageWindow(filename=os.path.basename(volfile),
                              image=self.projection,
                              dim=self.size, label="Projection")
        self.iw.updateData(flipud(self.projection.getData()))
        if self.showjPort:
            self.showjThread = Thread(target=self.listenShowJ)
            self.showjThread.daemon = True
            self.showjThread.start()
        self.iw.root.protocol("WM_DELETE_WINDOW", self.exitClient)
        self.iw.show()

    def rotate(self, rot, tilt, psi):
        self.fourierprojector.projectVolume(self.projection, rot, tilt, psi)
        self.projectionData = flipud(self.projection.getData())

        if hasattr(self, 'iw'):
            # sometimes is not created and rotate is called
            self.iw.updateData(self.projectionData)

    def exit(self):
        ChimeraClient.exit(self)
        if hasattr(self, "iw"):
            self.iw.root.destroy()

    def answer(self, msg):
        ChimeraClient.answer(self, msg)
        if msg == 'motion_stop':
            data = self.client.recv()  # wait for data
            printCmd('reading motion')
            self.motion = data
            printCmd('getting euler angles')
            rot, tilt, psi = xmipp.Euler_matrix2angles(self.motion)
            printCmd('calling rotate')
            self.rotate(rot, tilt, psi)

    def exitClient(self):  # close window before volume loaded
        if not self.listen:
            sys.exit(0)

    def initListenThread(self):
        self.listen_thread = Thread(target=self.listen)
        self.listen_thread.daemon = True
        self.listen_thread.start()

    def listenShowJ(self):
        self.serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.serversocket.bind(('', self.showjPort))
        # become a server socket
        self.serversocket.listen(1)

        while True:
            try:
                (clientsocket, address) = self.serversocket.accept()
                # should be a single message, so no loop
                msg = clientsocket.recv(1024)
                tokens = shlex.split(msg)
                cmd = tokens[0]
                if cmd == 'rotate':
                    rot = float(tokens[1])
                    tilt = float(tokens[2])
                    psi = float(tokens[3])

                    matrix = xmipp.Euler_angles2matrix(rot, tilt, psi)
                elif cmd == 'rotate_matrix':
                    matrixString = tokens[1]
                    matrix = ast.literal_eval(matrixString)
                    matrix = [matrix[0][:3], matrix[1][:3], matrix[2][:3]]
                self.client.send('rotate')
                self.client.send(matrix)
                clientsocket.close()

            except EOFError:
                print ('Lost connection to client')


def printCmd(cmd):
    pass
    # print cmd


def getVmdEnviron():
    """ Return the proper environ to launch VMD.
    VMD_HOME variable is read from the ~/.config/scipion.conf file.
    """
    environ = Environ(os.environ)
    environ.set('PATH', os.path.join(os.environ['VMD_HOME'], 'bin'),
                position=Environ.BEGIN)
    return environ


class VmdView(CommandView):
    """ View for calling an external command. """
    def __init__(self, vmdCommand, **kwargs):
        CommandView.__init__(self, 'vmd %s' % vmdCommand,
                             env=getVmdEnviron(), **kwargs)

    def show(self):
        runJob(None, '', self._cmd, env=getVmdEnviron())


class VmdViewer(Viewer):
    """ Wrapper to visualize PDB objects with VMD viewer. """
    _environments = [DESKTOP_TKINTER]
    # _targets = [PdbFile]

    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):
        cls = type(obj)

        if issubclass(cls, PdbFile):
            VmdView(obj.getFileName()).show()
            # FIXME: there is an asymetry between ProtocolViewer and Viewer.
            # For the first, the visualize method return a list of View's,
            # while for the second, the visualize method directly shows
            # the objects. (the first approach is preferable)
        else:
            raise Exception('VmdViewer.visualize: can not visualize class: %s'
                            % obj.getClassName())
