#!/usr/bin/python

import os

emPath = "http://scipionwiki.cnb.csic.es/files/scipion/software/em/"
brandeis = emPath + "brandeis.tgz"
eman = emPath + "eman2.pre2-1.linux64.tar.gz"
spider = emPath + "spiderweb.21.13.tar.gz"
relion = emPath + "relion-1.2.tar.bz2"

os.system("wget " + brandeis)
os.system("wget " + eman)
os.system("wget " + spider)
os.system("wget " + relion)