#!/usr/bin/env python
from platform import platform
from urllib import urlopen
#urlopen('http://localhost:8080/xmippWebServices/ReportInstallation?osystem='+platform())
#urlopen('http://localhost:8080/xmippWebServices/ReportInstallation?osystem=doNotStore')
#urlopen('http://glassfishdev.cnb.csic.es:8880/xmippWebServices/ReportInstallation?osystem='+platform())
#urlopen('http://glassfishprod.cnb.csic.es:8880/xmippWebServices/ReportInstallation?osystem='+platform())
a = urlopen('http://i2pc.cnb.csic.es/xmippWebServices/ReportInstallation?osystem='+platform())
#urlopen('http://i2pc.cnb.csic.es/xmippWebServices/ReportInstallation?osystem='+'kkkk')

