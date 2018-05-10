import requests
import os
import re
from importlib import import_module
import json
import pkg_resources
from pkg_resources import parse_version
from pyworkflow.utils.path import cleanPath
from pyworkflow import LAST_VERSION, OLD_VERSIONS
from install.funcs import Environment

REPOSITORY_URL = os.environ.get('SCIPION_PLUGIN_JSON', None) or os.environ['SCIPION_PLUGIN_REPO_URL']
PIP_BASE_URL = 'https://pypi.python.org/pypi'
PIP_CMD = 'python {}/pip install %(installSrc)s'.format(Environment.getPythonPackagesFolder())

versions = list(OLD_VERSIONS) + [LAST_VERSION]
if os.environ['SCIPION_SHORT_VERSION'] not in versions:
    SCIPION_VERSION = LAST_VERSION
else:
    SCIPION_VERSION = os.environ['SCIPION_SHORT_VERSION']

class PluginInfo(object):

    def __init__(self, pipName, name="", pluginSourceUrl="", remote=True, **kwargs):
        self.pipName = pipName
        self.name = name
        self.pluginSourceUrl = pluginSourceUrl
        self.remote = remote

        # things from pypi
        self.homePage = ""
        self.summary = ""
        self.author = ""
        self.email = ""
        self.compatibleReleases = {}
        self.latestRelease = ""

        # things we have when installed
        self.dirName = ""
        self.pipVersion = ""
        self.pipPath = ""
        self.emLink = ""
        self.binVersions = []
        self.pluginEnv = None

        if self.remote:
            self.setRemotePluginInfo()

        self.setLocalPluginInfo()  # get local info if installed

    def install(self):
        self.installPipModule()
        self.installBin()

    def hasPipPackage(self):
        try:
            pkg_resources.get_distribution(self.pipName)
            return True
        except Exception as e:
            return False

    def isInstalled(self):
        return self.hasPipPackage()

    def setRemotePluginInfo(self):
        reg = r'scipion-([\d.]*\d)'
        pipData = requests.get("%s/%s/json" % (PIP_BASE_URL, self.pipName))
        if pipData.ok:
            pipData = pipData.json()
            pluginInfo = pipData['info']
            releases = {}
            latestCompRelease = "0.0.0"
            for release, releaseData in pipData['releases'].iteritems():
                releaseData = releaseData[0]
                scipionVersions = [parse_version(v) for v in re.findall(reg, releaseData['comment_text'])]
                if len(scipionVersions) == 0:
                    print("WARNING: %s's release %s did not specify a compatible Scipion version" % (self.pipName,
                                                                                                     release))
                elif any([v <= parse_version(SCIPION_VERSION) for v in scipionVersions]):
                    if parse_version(latestCompRelease) < parse_version(release):
                        latestCompRelease = release
                    releases[release] = releaseData

            self.homePage = pluginInfo['home_page']
            self.summary = pluginInfo['summary']
            self.author = pluginInfo['author']
            self.email = pluginInfo['author_email']
            self.compatibleReleases = releases
            self.latestRelease = latestCompRelease

    def setLocalPluginInfo(self):
        if self.isInstalled():
            # metadata = json.loads(pkg_resources.get_distribution(self.pipName).get_metadata('metadata.json'))
            package = pkg_resources.get_distribution(self.pipName)
            keys = ['Name', 'Version', 'Summary', 'Home-page', 'Author', 'Author-email']
            pattern = r'(.*): (.*)'
            metadata = {}
            for line in package._get_metadata(package.PKG_INFO):
                match = re.match(pattern, line)
                if match:
                    key = match.group(1)
                    if key in keys:
                        metadata[key] = match.group(2)
                        keys.remove(key)
                        if not len(keys):
                            break

            self.pipVersion = metadata.get('Version', "")
            self.dirName = self.getDirName()
            self.pipPath = self.getPipPath()
            self.emLink = self.getEmPackagesLink()
            self.binVersions = self.getBinVersions()

            if not self.remote:
                # if we don't already have this info from remote, load it from metadata.json
                # self.homePage = metadata['extensions']['python.details'].get('project_urls', {}).get('Home', "")
                # self.summary = metadata['summary']
                # contacts = metadata['extensions']['python.details'].get('contacts', None)
                # self.author = contacts[0]['name'] if contacts else ""
                # self.email = contacts[0]['email'] if contacts else ""
                self.homePage = metadata.get('Home-page', "")
                self.summary = metadata.get('Summary', "")
                self.author = metadata.get('Author', "")
                self.email = metadata.get('Author-email', "")

    def getPluginObj(self):
        if os.path.exists(os.path.join(Environment.getPythonPackagesFolder(), self.dirName, 'plugin.py')):
            plugin = getattr(import_module('%s.plugin' % self.dirName), '_plugin')
        else:
            print("Warning: couldn't find _plugin in %s" % self.dirName)
            plugin = None
        return plugin

    def getInstallenv(self, envArgs=None):
        if envArgs is None:
            envArgs = []
        environment = Environment(args=envArgs)
        plugin = self.getPluginObj()
        if plugin:
            plugin.registerPluginBinaries(environment)
        return environment

    def getBinVersions(self):
        environment = self.getInstallenv()
        binVersions = [target.getName() for target in environment.getTargetList()]
        return binVersions

    def getDirName(self):
        # top level file is a file included in all pip packages that contains
        # the name of the package's top level directory
        return pkg_resources.get_distribution(self.pipName).get_metadata('top_level.txt').strip()

    def getPipPath(self):
        if self.dirName:
            return os.path.join(Environment.getPythonPackagesFolder(), self.dirName)
        else:
            return ""

    def getEmPackagesLink(self):
        if self.dirName:
            return os.path.join(Environment.getEmPackagesFolder(), self.dirName)
        else:
            return ""

    def installPipModule(self, version=""):
        environment = Environment()

        if not version:
            version = self.latestRelease
        elif version not in self.compatibleReleases:
            if self.compatibleReleases:
                print('%s version %s not compatible with current Scipion version %s.' % (self.pipName,
                                                                                         version,
                                                                                         SCIPION_VERSION))
                print("Please choose a compatible release: %s" % " ".join(self.compatibleReleases.keys()))

            else:
                print("%s has no compatible versions with current Scipion version %s." % (self.pipName,
                                                                                          SCIPION_VERSION))
            return False

        cmd = PIP_CMD % {'installSrc': self.pluginSourceUrl or "%s==%s" % (self.pipName, version)}

        pipModule = environment.addPipModule(self.pipName,
                                             target="%s*" % self.pipName.replace('-', '_'),
                                             pipCmd=cmd,
                                             ignoreDefaultDeps=True)
        environment.execute()

        # we already have a dir for the plugin:
        self.dirName = self.getDirName()

        # link pip module in em/packages
        if not os.path.exists(self.getEmPackagesLink()):
            os.symlink(self.getPipPath(), self.getEmPackagesLink())

        return True

    def installBin(self, args=None):
        environment = self.getInstallenv(envArgs=args)
        environment.execute()
        self.setLocalPluginInfo()

    def uninstallBins(self, binList=None):
        if binList is None:
            binList = self.binVersions

        binFolder = Environment.getEmFolder()
        for binVersion in binList:
            f = os.path.join(binFolder, binVersion)
            if os.path.exists(f):
                print('Removing %s binaries...' % binVersion)
                realPath = os.path.realpath(f)  # in case its a link
                cleanPath(f, realPath)
        return

    def uninstallPip(self):
        print('Removing %s plugin...' % self.pipName)
        try:
            from pip import main as pipmain
        except:
            from pip._internal import main as pipmain
        if os.path.exists(self.emLink):
            os.remove(self.emLink)
        pipmain(['uninstall', '-y', self.pipName])
        return

    def printBinInfo(self):
        env = self.getInstallenv()
        return env.printHelp()

class PluginRepository(object):

    def __init__(self, repoUrl=REPOSITORY_URL):
        self.repoUrl = repoUrl

    def getPlugins(self, pluginList=None, getPipData=False):

        pluginsJson = {}
        pluginDict = {}

        if os.path.isfile(self.repoUrl):
            with open(self.repoUrl) as f:
                pluginsJson = json.load(f)
        else:
            r = requests.get(self.repoUrl)
            if r.ok:
                pluginsJson = r.json()
            else:
                print("WARNING: Can't get scipion's plugin list, the plugin repository is not available")
                return pluginDict

        availablePlugins = pluginsJson.keys()

        if pluginList is None:
            targetPlugins = availablePlugins
        else:
            targetPlugins = set(availablePlugins).intersection(set(pluginList))
            if len(targetPlugins) < len(pluginList):
                wrongPluginNames = set(pluginList) - set(availablePlugins)
                print("WARNING - The following plugins didn't match available plugin names:")
                print(" ".join(wrongPluginNames))

        for pluginName in targetPlugins:
            pluginsJson[pluginName].update(remote=getPipData)
            pluginDict[pluginName] = PluginInfo(**pluginsJson[pluginName])
        return pluginDict


    def printPluginInfo(self, withBins=False):
        printStr = ""
        pluginDict = self.getPlugins()
        if pluginDict:
            withBinsStr = "Installed plugins and their binaries" if withBins else "Available plugins"
            printStr += ("%s: "
                         "([ ] not installed, [X] seems already installed)\n\n" % withBinsStr)
            keys = sorted(pluginDict.keys())
            for name in keys:
                plugin = pluginDict[name]
                if withBins and not plugin.isInstalled():
                    continue
                printStr += "%23s " % name
                vInfo = '[%s]' % ('X' if plugin.isInstalled() else ' ')
                printStr += '%13s' % vInfo
                printStr += "\n"
                if withBins:
                    printStr += plugin.printBinInfo().split('\n', 1)[1]
        else:
            printStr = "List of available plugins in plugin repository unaccessible at this time."
        return printStr
