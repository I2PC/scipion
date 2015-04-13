"""
    http://www.scons.org/wiki/UnTarBuilder

    env.UnTar(source = 'apr-1.3.3')

"""

import os
from SCons.Builder import Builder
from SCons.Action import Action
#from SCons.Node.FS import File


def tarContentsEmitter(target, source, env):
    import tarfile
    sourceTar = tarfile.open(source[0].name,'r')
    tarContents = sourceTar.getmembers()
    tarFileContents = filter(lambda tarEntry: tarEntry.isfile(), tarContents)
    newTargets = map(tarInfoToNode, tarFileContents)
    sourceTar.close()
    return (newTargets, source)

def tarInfoToNode(tarInfoObject):
    return tarInfoObject.name
    #return File(tarInfoObject.name)

def UnTar(target, source, env):
    # Code to build "target" from "source" here
    import tarfile
    sourceTar = tarfile.open(source[0].name,'r')
    sourceTar.extractall()
    sourceTar.close()
    return None

def UnTarString(target, source, env):
    """ Information string for UnTar """
    return 'Extracting %s' % os.path.basename (str (source[0]))

unTarBuilder = Builder(action=Action(UnTar, UnTarString),
                       src_suffix='.tar.gz',
                       emitter=tarContentsEmitter)

def generate(env):
    env.Append(BUILDERS = {'UnTar' : unTarBuilder })

def exists(env):
    return True

