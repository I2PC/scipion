# Run django site in apache with wsgi
import os, sys
sys.path.append('/home/scipion/pyworkflow-code/pyworkflow/web')
sys.path.append('/home/scipion/pyworkflow-code')

# Some important environment variables, like LD_LIBRARY_PATH, must be set before running this script
# They are set in /etc/apache2/envvars

# Where to look for xmipp python code
xmipp_home='/home/scipion/xmipp'
sys.path.append(xmipp_home + '/lib')
sys.path.append(xmipp_home + '/protocols')
sys.path.append(xmipp_home + '/applications/tests/pythonlib')
sys.path.append(xmipp_home + '/lib/python2.7/site-packages')

os.environ['DJANGO_SETTINGS_MODULE'] = 'pages.settings'

import django.core.handlers.wsgi

application = django.core.handlers.wsgi.WSGIHandler()
