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

    ####################### Install funcs ############################

    def install(self):
        """Installs both pip module and default binaries of
        the plugin"""
        self.installPipModule()
        self.installBin()

    def hasPipPackage(self):
        """Checks if the current plugin is installed via pip"""
        try:
            pkg_resources.get_distribution(self.pipName)
            return True
        except Exception as e:
            return False

    def isInstalled(self):
        """Checks if the current plugin is installed (i.e. has pip package).
        NOTE: we might wanna change definition of isInstalled, hence the extra function."""
        return self.hasPipPackage()

    def installPipModule(self, version=""):
        """Installs the version specified of the pip plugin, as long as it is compatible
        with the current Scipion version. If no version specified, will install latest
        compatible one."""
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
        reloadPkgRes = self.isInstalled()  # check if we're doing a version change of an already installed plugin
        environment.execute()
        if reloadPkgRes:
            # if plugin was already installed, pkg_resources has the old one so it needs a reload
            reload(pkg_resources)
        # we already have a dir for the plugin:
        self.dirName = self.getDirName()

        # link pip module in em/packages
        try:
            os.symlink(self.getPipPath(), self.getEmPackagesLink())
        except OSError:
            if os.path.islink(self.getEmPackagesLink()):
                return True
            else:
                raise

        return True

    def installBin(self, args=None):
        """Install binaries of the plugin. Args is the list of args to be passed to
        the install environment."""
        environment = self.getInstallenv(envArgs=args)
        environment.execute()
        self.setLocalPluginInfo()

    def uninstallBins(self, binList=None):
        """Install binaries of the plugin.
        - binList: if  given, will install the binaries in it. The binList
                   may contain strings with only the name of the binary or
                   name and version in the format name-version"""
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
        """Removes pip package from site-packages and removes link in
        pyworkflow/em/packages"""
        print('Removing %s plugin...' % self.pipName)
        try:
            from pip import main as pipmain
        except:
            from pip._internal import main as pipmain
        if os.path.exists(self.emLink):
            os.remove(self.emLink)
        pipmain(['uninstall', '-y', self.pipName])
        return

    ####################### Remote data funcs ############################

    def getPipJsonData(self):
        """"Request json data from pypi, return json content"""
        pipData = requests.get("%s/%s/json" % (PIP_BASE_URL, self.pipName))
        if pipData.ok:
            pipData = pipData.json()
            return pipData
        else:
            print("Warning: Couldn't get remote plugin data for %s" % self.pipName)
            return {}

    def getCompatiblePipReleases(self, pipJsonData=None):
        """Get pip releases of this plugin that are compatible with
         current Scipion version. Returns dict with all compatible releases and
         a special key "latest" with the most recent one."""

        if pipJsonData is None:
            pipJsonData = self.getPipJsonData()

        reg = r'scipion-([\d.]*\d)'
        releases = {}
        latestCompRelease = "0.0.0"

        for release, releaseData in pipJsonData['releases'].iteritems():
            releaseData = releaseData[0]
            scipionVersions = [parse_version(v) for v in re.findall(reg, releaseData['comment_text'])]
            if len(scipionVersions) == 0:
                print("WARNING: %s's release %s did not specify a compatible Scipion version" % (self.pipName,
                                                                                                 release))
            elif any([v <= parse_version(SCIPION_VERSION.strip('v')) for v in scipionVersions]):
                if parse_version(latestCompRelease) < parse_version(release):
                    latestCompRelease = release
                releases[release] = releaseData

        releases['latest'] = latestCompRelease

        return releases

    def setRemotePluginInfo(self):
        """Sets value for the attributes that need to be obtained from pypi"""
        reg = r'scipion-([\d.]*\d)'
        pipData = self.getPipJsonData()
        if not pipData:
            return
        info = pipData['info']
        releases = self.getCompatiblePipReleases(pipJsonData=pipData)

        self.homePage = info['home_page']
        self.summary = info['summary']
        self.author = info['author']
        self.email = info['author_email']
        self.compatibleReleases = releases
        self.latestRelease = releases['latest']

    ####################### Local data funcs ############################

    def setLocalPluginInfo(self):
        """Sets value for the attributes that can be obtained locally if the
        plugin is installed."""
        if self.isInstalled():
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

            if not self.remote:  # only do this if we don't already have it from remote
                self.homePage = metadata.get('Home-page', "")
                self.summary = metadata.get('Summary', "")
                self.author = metadata.get('Author', "")
                self.email = metadata.get('Author-email', "")

    def getPluginObj(self):
        """ Tries to find the _plugin object in plugin.py file."""
        if os.path.exists(os.path.join(Environment.getPythonPackagesFolder(), self.dirName, 'plugin.py')):
            plugin = getattr(import_module('%s.plugin' % self.dirName), '_plugin')
        else:
            print("Warning: couldn't find _plugin in %s" % self.dirName)
            plugin = None
        return plugin

    def getInstallenv(self, envArgs=None):
        """Reads the registerPluginBinaries function from plugin.py and returns an
        Environment object with the plugin's binaries."""
        if envArgs is None:
            envArgs = []
        environment = Environment(args=envArgs)
        plugin = self.getPluginObj()
        if plugin:
            plugin.registerPluginBinaries(environment)
        return environment

    def getBinVersions(self):
        """Get list with names of binaries of this plugin"""
        environment = self.getInstallenv()
        binVersions = [target.getName() for target in environment.getTargetList()]
        return binVersions

    def getDirName(self):
        """Get the name of the folder that contains the plugin code itself (e.g. to import
        the _plugin object.)"""
        # top level file is a file included in all pip packages that contains
        # the name of the package's top level directory
        return pkg_resources.get_distribution(self.pipName).get_metadata('top_level.txt').strip()

    def getPipPath(self):
        """Get path of the plugin in site packages folder"""
        if self.dirName:
            return os.path.join(Environment.getPythonPackagesFolder(), self.dirName)
        else:
            return ""

    def getEmPackagesLink(self):
        """Get path of the plugin link in Em packages folder"""
        if self.dirName:
            return os.path.join(Environment.getEmPackagesFolder(), self.dirName)
        else:
            return ""

    def printBinInfoStr(self):
        """Returns string with info of binaries installed to print in console
        with flag --help"""
        env = self.getInstallenv()
        return env.printHelp()

class PluginRepository(object):

    def __init__(self, repoUrl=REPOSITORY_URL):
        self.repoUrl = repoUrl

    def getPlugins(self, pluginList=None, getPipData=False):
        """Reads available plugins from self.repoUrl and returns a dict with
        PluginInfo objects. Params:
        - pluginList: A list with specific plugin pip-names we want to get.
        - getPipData: If true, each PluginInfo object will try to get the data
        of the plugin from pypi."""

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
                print("WARNING: Can't get Scipion's plugin list, the plugin repository is not available")
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

    def printPluginInfoStr(self, withBins=False, withUpdates=False):
        """Returns string to print in console which plugins are installed.
        - withBins: If true, will add binary info for the plugins installed
        - with Updates: If true, will check if the installed plugins have new
                    releases."""
        def ansi(n):
            """Return function that escapes text with ANSI color n."""
            return lambda txt: '\x1b[%dm%s\x1b[0m' % (n, txt)

        black, red, green, yellow, blue, magenta, cyan, white = map(ansi, range(30, 38))

        printStr = ""
        pluginDict = self.getPlugins(getPipData=withUpdates)
        if pluginDict:
            withBinsStr = "Installed plugins and their binaries" if withBins else "Available plugins"
            printStr += ("%s: "
                         "([ ] not installed, [X] seems already installed)\n\n" % withBinsStr)
            keys = sorted(pluginDict.keys())
            for name in keys:
                plugin = pluginDict[name]
                if withBins and not plugin.isInstalled():
                    continue
                printStr += "%23s" % name
                vInfo = '%s [%s]' % (plugin.pipVersion, 'X' if plugin.isInstalled() else ' ')
                printStr += '%15s' % vInfo
                if withUpdates and plugin.isInstalled():
                    if plugin.latestRelease != plugin.pipVersion:
                        printStr += yellow('\t(%s available)' % plugin.latestRelease)
                printStr += "\n"
                if withBins:
                    printStr += plugin.printBinInfoStr().split('\n', 1)[1]
        else:
            printStr = "List of available plugins in plugin repository inaccessible at this time."
        return printStr
