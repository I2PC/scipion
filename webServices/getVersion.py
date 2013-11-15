#!/usr/bin/env python
from urllib import urlopen

#urlopen('http://localhost:8080/xmippWebServices/GetVersion')
#urlopen('http://glassfishdev.cnb.csic.es:8880/xmippWebServices/GetVersion')
#urlopen('http://glassfishprod.cnb.csic.es:8880/xmippWebServices/GetVersion')
#urlopen('http://i2pc.cnb.csic.es/xmippWebServices/GetVersion')
a = urlopen('http://localhost:8080/xmippWebServices/GetVersion')
exec (a.read())
print git_date, version, git_hash

