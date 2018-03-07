import requests
import os
from glob import glob
from pyworkflow.utils.path import cleanPath
from install.funcs import Environment, Command, Link
import pip

REPOSITORY_URL = 'http://localhost:8000/plugins.json'
PIP_CMD = 'python {}/pip install %(installSrc)s'.format(Environment.getPythonPackagesFolder())

class PluginInfo(object):

    def __init__(self, name, dirName,  pipName="", pluginVersions=None,
                 binaryVersions=None, pluginSourceUrl=""):
        self.name = name
        self.dirName = dirName
        self.pluginVersions = pluginVersions
        self.pluginSourceUrl = pluginSourceUrl
        self.pipName = pipName
        self.binaryVersions = binaryVersions
        self.pluginPipPath = os.path.join(Environment.getPythonPackagesFolder(), self.dirName)
        self.pluginEMLink = os.path.join(Environment.getEmPackagesFolder(), self.dirName)


    def install(self):
        self.installPipModule()
        self.installBin()

    def installPipModule(self):
        environment = Environment()
        cmd = PIP_CMD % {'installSrc': self.pluginSourceUrl}

        pipModule = environment.addPipModule(self.dirName,
                                             pipCmd=cmd,
                                             ignoreDefaultDeps=True,)

        # link pip module to em/packages
        pipModule.addCommand(Command(environment, Link(self.pluginEMLink, self.pluginPipPath),
                                     targets=[self.pluginEMLink]))
        environment.execute()

    def installBin(self, **kwargs):
        environment = Environment(**kwargs)
        # install binaries
        from importlib import import_module
        if os.path.exists(os.path.join(environment.getPythonPackagesFolder(), self.dirName, 'plugin.py')):
            plugin = getattr(import_module('%s.plugin' % self.dirName), '_plugin')
            target = plugin.registerFunction(environment)
        else:
            print("Cant find plugin.py for %s" % self.dirName)

        environment.execute()

    def uninstallBin(self, version=""):

        if version:
            name = Environment._getExtName(self.name, version)
        else:
            name = self.name
        files = glob("%s*" % os.path.join(Environment.getEmFolder(), name))
        if files:
            print('Removing %s binaries...' % name)
            for f in files:
                realPath = os.path.realpath(f)  # in case its a link
                cleanPath(f, realPath)

        return

    def uninstallPip(self):
        print('Removing %s plugin...' % self.pipName)
        if os.path.exists(self.pluginEMLink):
            os.remove(self.pluginEMLink)
        pip.main(['uninstall', '-y', self.pipName])

        return

class PluginRepository(object):

    def __init__(self, environment=None, repoUrl=REPOSITORY_URL):
        self.repoUrl = repoUrl
        self.pluginDict = self.getAvailablePlugins()
        self.environment = environment or Environment()

    def getAvailablePlugins(self):
        r = requests.get(self.repoUrl)
        if r.ok:
            pluginDict = {k: PluginInfo(**v) for k, v in r.json().iteritems()}
            return pluginDict
        else:
            print("WARNING: Can't get scipion's plugin list, the plugin repository isn't unavailable")
            return {}

    def printPluginInfo(self):
        printStr = ""
        if self.pluginDict:
            printStr += ("Available plugins: "
                         "([ ] not installed, [X] seems already installed)\n\n")
            keys = sorted(self.pluginDict.keys())
            for name in keys:
                pVersions = getattr(self.pluginDict[name], 'binaryVersions')
                printStr += "%15s " % name
                for version in pVersions:
                    installed = self.environment._isInstalled(name, version)
                    vInfo = '%s [%s]' % (version, 'X' if installed else ' ')
                    printStr += '%13s' % vInfo

                printStr += "\n"
        else:
            printStr = "List of available plugins in plugin repository unaccessible at this time."
        return printStr



