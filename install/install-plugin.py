import os
import sys
from install.funcs import Environment
from pyworkflow.em import PACKAGES_PATH
#  ************************************************************************
#  *                                                                      *
#  *                       External (EM) Plugins                          *
#  *                                                                      *
#  ************************************************************************

args = sys.argv
env = Environment(args=sys.argv)

for name, plugin in env.getPlugins(PACKAGES_PATH).iteritems():
    plugin.registerFunction(env)

env.execute()
