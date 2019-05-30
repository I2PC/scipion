# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import os

try:  # python 2
    import Tkinter as tk
    import tkFont
    import ttk
except ImportError:  # Python 3
    import tkinter as tk
    import tkinter.font as tkFont
    import tkinter.ttk as ttk

from pyworkflow.em.convert import ImageHandler
from pyworkflow.viewer import View, Viewer, CommandView, DESKTOP_TKINTER, ProtocolViewer
from pyworkflow.utils import Environ, runJob, importFromPlugin, getFreePort

import showj


class DataView(View):
    """ Wrapper the arguments to showj (either web or desktop).
    Also useful to visualize images that are not objects,
    e.g.: dark or gain images"""
    def __init__(self, path, viewParams={}, **kwargs):
        View.__init__(self)
        self._memory = showj.getJvmMaxMemory()
        self._loadPath(path)
        self._env = kwargs.get('env', {})
        self._viewParams = viewParams

    def setMemory(self, memory):
        self._memory = memory

    def getViewParams(self):
        """ Give access to the viewParams dict. """
        return self._viewParams

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
        showj.runJavaIJapp(self._memory,
                           'xmipp.viewer.scipion.ScipionViewer',
                           self.getShowJParams(), env=self._env)

    def getShowJParams(self):
        tableName = '%s@' % self._tableName if self._tableName else ''
        params = '-i "%s%s"' % (tableName, self._path)
        for key, value in self._viewParams.items():
            params = "%s --%s %s" % (params, key, value)

        return params

    def getShowJWebParams(self):
    
        parameters = {
            showj.MODE,  # FOR MODE TABLE OR GALLERY
            showj.VISIBLE,
            showj.ZOOM,
            showj.ORDER,
            showj.RENDER,
            showj.SORT_BY
        }
        
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

        def existingLabels(labelList):

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
                  '_xmipp_ctfmodel_halfplane', '_micObj.plotGlobal._filename'
                 ]
    EXTRA_LABELS = ['_ctffind4_ctfResolution', '_gctf_ctfResolution',
                    '_ctffind4_ctfPhaseShift', '_gctf_ctfPhaseShift',
                    '_ctftilt_tiltAxis', '_ctftilt_tiltAngle',
                    '_xmipp_ctfCritFirstZero',
                    '_xmipp_ctfCritCorr13', '_xmipp_ctfCritIceness','_xmipp_ctfCritFitting',
                    '_xmipp_ctfCritNonAstigmaticValidty',
                    '_xmipp_ctfCritCtfMargin', '_xmipp_ctfCritMaxFreq',
                    '_xmipp_ctfCritPsdCorr90', '_xmipp_ctfVPPphaseshift'
                   ]

    def __init__(self, project, ctfSet, other='', **kwargs):
        first = ctfSet.getFirstItem()

        def existingLabels(labelList):
            return ' '.join([l for l in labelList if first.hasAttributeExt(l)])

        psdLabels = existingLabels(self.PSD_LABELS)
        extraLabels = existingLabels(self.EXTRA_LABELS)
        labels =  'id enabled %s _defocusU _defocusV ' % psdLabels
        labels += '_defocusAngle _defocusRatio '
        labels += '_phaseShift _resolution _fitQuality %s ' % extraLabels
        labels += ' _micObj._filename'

        viewParams = {showj.MODE: showj.MODE_MD,
                      showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.ZOOM: 50
                     }

        if psdLabels:
            viewParams[showj.RENDER] = psdLabels

        if ctfSet.isStreamOpen():
            viewParams['dont_recalc_ctf'] = ''

        def _anyAttrStartsBy(obj, prefix):
            """ Return True if any of the attributes of this object starts
            by the provided prefix.
            """
            return any(attrName.startswith(prefix)
                       for attrName, _ in obj.getAttributesToStore())

        if _anyAttrStartsBy(first, '_ctffind4_ctfResolution'):
            gviewer = importFromPlugin('grigoriefflab.viewers', '')
            viewParams[showj.OBJCMDS] = "'%s'" % gviewer.OBJCMD_CTFFIND4

        elif _anyAttrStartsBy(first, '_gctf'):
            OBJCMD_GCTF = importFromPlugin('gctf.viewers', 'OBJCMD_GCTF')
            viewParams[showj.OBJCMDS] = "'%s'" % OBJCMD_GCTF

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
    MODE_AUTOMATIC = 'Automatic'

    def __init__(self, project, path, outputdir, protocol, pickerProps=None,
                 inTmpFolder=False, **kwargs):
        DataView.__init__(self, path, **kwargs)
        self.project = project
        self.outputdir = outputdir
        self.protocol = protocol
        self.pickerProps = pickerProps
        self.inTmpFolder = inTmpFolder
        self.mode = kwargs.get('mode', None)

    def show(self):
        return showj.launchSupervisedPickerGUI(self._path, self.outputdir,
                                               self.protocol, mode=self.mode,
                                               pickerProps=self.pickerProps,
                                               inTmpFolder=self.inTmpFolder)


class ImageView(View):
    """ Customized ObjectView for SetOfClasses. """
    def __init__(self, imagePath, **kwargs):
        View.__init__(self)
        self._imagePath = os.path.abspath(imagePath)

    def getImagePath(self):
        return self._imagePath


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
                 fontSize=16, padding=10,
                 fontFamily='monospace'):
        # get new widget that has as parent the top level window and set title
        win = tk.Toplevel()
        if title:
            win.wm_title(title)

        # frame to place all other widgets
        frame = tk.Frame(win)

        # make font a little bigger
        # TODO: font size should be general
        font = tkFont.Font(family=fontFamily, size=fontSize)
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
