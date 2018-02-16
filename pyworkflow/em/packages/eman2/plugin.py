from pyworkflow.plugin import Plugin
from pyworkflow.em.packages.eman2 import EMAN_DIR_VAR
import os

VARS = {EMAN_DIR_VAR: os.path.join(os.environ['SCIPION_SOFTWARE'], 'eman-2.12')}


def registerPluginBinaries(env):
    eman2_commands = [('./eman2-installer',
                       'eman2.*rc')]

    env.addPackage('eman', version='2.11',
                   url='http://scipion.cnb.csic.es/downloads/scipion/software/em/eman2.11.linux64.tgz',
                   tar='eman2.11.linux64.tgz',
                   commands=eman2_commands)

    env.addPackage('eman', version='2.12',
                   url='http://scipion.cnb.csic.es/downloads/scipion/software/em/eman2.12.linux64.tgz',
                   tar='eman2.12.linux64.tgz',
                   commands=eman2_commands)


_plugin = Plugin('eman2',
                 version=2,
                 configVars=VARS,
                 logo="eman2_logo2.png",
                 references=['Tang2007'],
                 registerFunction=registerPluginBinaries)

