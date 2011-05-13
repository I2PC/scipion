XMIPP_BASE='/home/xmipp/XMIPP'
XMIPP_HOME=XMIPP_BASE+'/xmipp'
XMIPP_TEST=XMIPP_HOME + '/applications/tests'
XMIPP_LOGS=XMIPP_BASE + '/Logs'
SERVER = 'www.cnb.csic.es'

fromaddr = "xmipp@bioweb.cnb.csic.es"
toaddrs  = ["roberto.marabini.cnb@gmail.com","jcuenca@cnb.csic.es","delarosatrevin@gmail.com","coss@cnb.csic.es","jvega@cnb.csic.es","joton@cnb.csic.es","alejandro.e.rey@gmail.com"]
#toaddrs  = ["roberto.marabini.cnb@gmail.com"]
subject = "xmipp Compilation"
testNames=['test_fftw'\
           'test_filters'\
           'xmipp_test_funcs',\
	   'xmipp_test_image',\
	   'xmipp_test_image_generic',\
	   'xmipp_test_metadata',\
	   'xmipp_test_multidim']
