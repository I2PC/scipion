# AutoConfig Builder:  Runs ./configure inside a directory.
#
# Parameters:
#    AutoConfigParams -- Sequence of parameter strings to include on the
#                        configure command line.
#                        Default: []
#    AutoConfigTarget -- File that configure will create.
#                        Default: "config.h"
#    AutoConfigSource -- File that configure depends on.
#                        Default: "Makefile.in"
#    AutoConfigStdOut -- File where the output will be written to.
#                        Default: None

from os.path import join, dirname
import subprocess

from SCons.Script import Exit


def parms(target, source, env):
    """Assemble various AutoConfig parameters."""

    workdir = dirname(str(source[0]))
    params = env.get('AutoConfigParams', [])
    if not isinstance(params, list):
        print 'AutoConfigParams must be a sequence'
        Exit(1)
    targetfile = env.get('AutoConfigTarget', 'config.h')
    sourcefile = env.get('AutoConfigSource', 'Makefile.in')
    out = env.get('AutoConfigStdOut')
    return (workdir, params, targetfile, sourcefile, out)


def message(target, source, env):
    """Return a pretty AutoConfig message."""

    dirx, params, targetfile, sourcefile, out = parms(target, source, env)

    if 'AUTOCONFIGCOMSTR' in env:
        msg = env.subst(env['AUTOCONFIGCOMSTR'],
                        target=target, source=source, raw=1)
        return '%s > %s' % (msg, out)

    return 'cd %s && ./configure %s' % (dirx, ' '.join(params))


def emitter(target, source, env):
    """Remap the source & target to path/$AutoConfigSource and path/$AutoConfigTarget."""

    dirx, params, targetfile, sourcefile, out = parms(target, source, env)

    # NOTE: Using source[0] instead of target[0] for the target's path!
    # If there's only one . in the source[0] value, then Scons strips off the
    # extension when it determines the target[0] value.  For example,
    #    AutoConfig('foo.blah')
    # sets
    #    source[0] = 'foo.blah'
    #    target[0] = 'foo'
    # (SCons does NOT do this if source[0] has more than one . )
    # Since this emitter expects the incoming source[0] value to be a directory
    # name, we can use it here for the rewritten target[0].

    return ([ join(str(source[0]), targetfile) ],
            [ join(str(source[0]), sourcefile) ])


def builder(target, source, env):
    """Run ./configure in a directory."""

    dirx, params, targetfile, sourcefile, out = parms(target, source, env)

    if 'AUTOCONFIGCOMSTR' in env and out is not None:
        fout = open(out, 'w+')
    else:
        fout = None

    return subprocess.call(['./configure'] + params, cwd=dirx,
                           stdout=fout, stderr=fout)


def generate(env, **kwargs):
    env['BUILDERS']['AutoConfig'] = env.Builder(
        action=env.Action(builder, message),
        emitter=emitter, single_source=True)


def exists(env):
    return True
