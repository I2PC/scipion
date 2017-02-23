#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import time
import psutil
import sys

# By default check xmipp programs
prefix = sys.argv[1] if len(sys.argv) > 1 else 'xmipp'
wait = int(sys.argv[2]) if len(sys.argv) > 2 else 24 * 60 * 60 # one day

print "Check program starting with: %s" % prefix
print "Check during %s seconds" % wait

for i in range(wait):
    c = 0
    for proc in psutil.process_iter():
        if proc.name().startswith(prefix):
            c += 1
    print "count: ", c
    time.sleep(1)

