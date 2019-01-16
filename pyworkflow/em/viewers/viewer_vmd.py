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

import pyworkflow.utils as pwutils
import pyworkflow.viewer as pwviewer
from pyworkflow.em.data import AtomStruct


class Vmd:
    """ Help class to run VMD and manage its environment. """

    @classmethod
    def getEnviron(cls):
        """ Return the proper environ to launch VMD.
        VMD_HOME variable is read from the ~/.config/scipion.conf file.
        """
        environ = pwutils.Environ(os.environ)
        environ.set('PATH', os.path.join(os.environ['VMD_HOME'], 'bin'),
                    position=pwutils.Environ.BEGIN)
        return environ


class VmdView(pwviewer.CommandView):
    """ View for calling an external command. """
    def __init__(self, vmdCommand, **kwargs):
        pwviewer.CommandView.__init__(self, 'vmd %s' % vmdCommand,
                             env=Vmd.getEnviron(), **kwargs)

    def show(self):
        pwutils.runJob(None, '', self._cmd, env=Vmd.getEnviron())


class VmdViewer(pwviewer.Viewer):
    """ Wrapper to visualize PDB objects with VMD viewer. """
    _environments = [pwviewer.DESKTOP_TKINTER]
    # _targets = [AtomStruct]

    def __init__(self, **args):
        pwviewer.Viewer.__init__(self, **args)

    def visualize(self, obj, **args):
        cls = type(obj)

        if issubclass(cls, AtomStruct):
            VmdView(obj.getFileName()).show()
            # FIXME: there is an asymmetry between ProtocolViewer and Viewer.
            # For the first, the visualize method return a list of View's,
            # while for the second, the visualize method directly shows
            # the objects. (the first approach is preferable)
        else:
            raise Exception('VmdViewer.visualize: can not visualize class: %s'
                            % obj.getClassName())
