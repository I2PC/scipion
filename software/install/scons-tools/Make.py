# Make Builder: Runs make.
#
# Parameters:
#    MakePath -- SCons Dir node representing the directory in which to run make.  REQUIRED.
#    MakeCmd -- The 'make' executable to run.
#               Default: make
#    MakeEnv -- Dictionary of variables to set in the make execution environment.
#               Default: none
#    MakeOpts -- Options to pass on the make command line.
#                Default: none
#    MakeOneThread -- Don't pass any -j option to make.
#                     Default: False
#    MakeTargets -- String of space-seperated targets to pass to make
#                   Default: ""

import os
import subprocess

from SCons.Script import Exit, GetOption


def parms(target, source, env):
    """Assemble various Make parameters."""

    if 'MakePath' not in env:
        print "Make builder requires MakePath variable"
        Exit(1)

    make_path = env.subst(str(env['MakePath']))

    make_cmd = 'make'
    if 'MakeCmd' in env:
        make_cmd = env.subst(env['MakeCmd'])
    elif 'MAKE' in env:
        make_cmd = env.subst(env['MAKE'])

    make_env = None
    if env.get('CROSS_BUILD'):
        make_env = env['CROSS_ENV']
    if 'MakeEnv' in env:
        if make_env == None:
            make_env = {}
        else:
            # We're appending to an existing dictionary, so create a copy
            # instead of appending to the original env['CROSS_ENV']
            make_env = env['CROSS_ENV'][:]
        for (k,v) in env['MakeEnv'].items():
            make_env[k] = v

    make_opts = None
    if 'MakeOpts' in env:
        make_opts = env.subst(env['MakeOpts'])

    make_jobs = GetOption('num_jobs')
    if 'MakeOneThread' in env and env['MakeOneThread']:
        make_jobs = 1

    make_targets = None
    if 'MakeTargets' in env:
        make_targets = env.subst(env['MakeTargets'])

    out = env.get('MakeStdOut')

    return (make_path, make_env, make_targets, make_cmd, make_jobs, make_opts, out)


def message(target, source, env):
    """Return a pretty Make message"""

    (make_path,
     make_env,
     make_targets,
     make_cmd,
     make_jobs,
     make_opts,
     out) = parms(target, source, env)

    myenv = env.Clone()
    # Want to use MakeTargets in the MAKECOMSTR, but make it pretty first.
    if 'MakeTargets' in myenv:
        myenv['MakeTargets'] += ' '
    else:
        myenv['MakeTargets'] = ''

    if 'MAKECOMSTR' in myenv:
        return myenv.subst(myenv['MAKECOMSTR'],
                           target=target, source=source, raw=1) + " > %s " % out

    msg = 'cd ' + make_path + ' &&'
    if make_env != None:
        for k, v in make_env.iteritems():
            msg += ' ' + k + '=' + v
    msg += ' ' + make_cmd
    if make_jobs > 1:
        msg += ' -j %d' % make_jobs
    if make_opts != None:
        msg += ' ' + ' '.join(make_opts)
    if make_targets != None:
        msg += ' ' + make_targets
    return msg


def builder(target, source, env):
    """Run make in a directory."""

    (make_path,
     make_env,
     make_targets,
     make_cmd,
     make_jobs,
     make_opts,
     out) = parms(target, source, env)

    # Make sure there's a directory to run make in
    if len(make_path) == 0:
        print 'No path specified'
    if not os.path.exists(make_path):
        print 'Path %s not found' % make_path

    # Build up the command and its arguments in a list
    fullcmd = [ make_cmd ]

    if make_jobs > 1:
        fullcmd += [ '-j', str(make_jobs) ]

    if make_opts:
        fullcmd += make_opts

    if make_targets:
        fullcmd += make_targets.split()

    # Capture the make command's output, unless we're verbose
    if out is not None:
        fout = open(out, 'w+')
    else:
        fout = None

    # Make!
    make = subprocess.Popen(fullcmd, cwd=make_path,
                            stdout=fout, stderr=fout,
                            env=make_env)

    # Some subprocesses don't terminate unless we communicate with them
    output = make.communicate()[0]
    return make.returncode


def generate(env, **kwargs):
    env['BUILDERS']['Make'] = env.Builder(action=env.Action(builder, message))


def exists(env):
    if env.WhereIs(env.subst('$MAKE')) != None:
        return True
    return False
