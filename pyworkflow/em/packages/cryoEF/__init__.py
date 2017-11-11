# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

from bibtex import _bibtex # Load bibtex dict with references
from convert import getEnviron

_logo = "cryoEF_logo.png"
_references = ['naydenova2017']

from protocol_cryoef import ProtCryoEF
from convert import *
from viewer_cryoef import CryoEFViewer
from wizard import cryoEFMaskDiameterWizard


_environ = getEnviron()


def validateInstallation():
    """ This function will be used to check if package is properly installed. """
    missingPaths = ["%s: %s" % (var, _environ[var])
                    for var in [CRYOEF_HOME_VAR]
                    if not os.path.exists(_environ[var])]

    if missingPaths:
        return ["Missing variables:"] + missingPaths
    else:
        return []  # No errors
