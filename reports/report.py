#!/usr/bin/env python
from platform import platform
from urllib import urlopen
#urlopen('http://localhost:8080/WebApplication3/NewServlet?osystem='+platform())
urlopen('http://localhost:8080/WebApplication3/NewServlet?osystem=doNotStore')
#urlopen('http://glassfishdev.cnb.csic.es:8880/WebApplication3/NewServlet?osystem='+platform())
#urlopen('http://glassfishprod.cnb.csic.es:8880/WebApplication3/NewServlet?osystem='+platform())
#urlopen('http://i2pc.cnb.csic.es/xmipp/NewServlet?osystem='+platform())

