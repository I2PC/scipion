# Run django site in apache with wsgi
# Migrated to scipion script
import os, sys
import syslog
# sys.path.append('/home/scipionweb/pyworkflow-code/pyworkflow/web')
# sys.path.append('/home/scipionweb/pyworkflow-code')

# Some important environment variables, like LD_LIBRARY_PATH, must be set before running this script
# They are set in /etc/apache2/envvars

# Where to look for xmipp python code
# xmipp_home='/home/scipionweb/pyworkflow-code/software/em/xmipp'
# sys.path.append(xmipp_home + '/lib')
# sys.path.append(xmipp_home + '/protocols')
# sys.path.append(xmipp_home + '/applications/tests/pythonlib')
# sys.path.append(xmipp_home + '/lib/python2.7/site-packages')


# SCIPION_HOME = os.environ.get('SCIPION_HOME', os.path.join('home','scipionweb', 'pyworkflow'))

# os.environ.update(module_globals["VARS"])


# syslog.syslog(os.environ.get('LD_LIBRARY_PATH'))

os.environ['DJANGO_SETTINGS_MODULE'] = 'pages.settings'
import django.core.handlers.wsgi

application = django.core.handlers.wsgi.WSGIHandler()
