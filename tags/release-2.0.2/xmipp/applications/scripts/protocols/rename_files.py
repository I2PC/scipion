#!/usr/bin/env python
import fileinput, glob, string, sys, os
from os.path import join
# replace a string in multiple files
#filesearch.py


def searchreplace(path,search,replace,exts=None):

	import fileinput, glob, string, sys, os
	from os.path import join
	# replace a string in multiple files
	#filesearch.py
	files = glob.glob(path + "/*")
	if files is not []:
		for file in files:
			if os.path.isfile(file):
				if exts is None or exts.count(os.path.splitext(file)[1]) is not 0:
					for line in fileinput.input(file,inplace=1):
						lineno = 0
						lineno = string.find(line, search)
						#sys.stderr.write('\n'+line)
						#sys.stderr.write(search)
						#sys.stderr.write(str(lineno)+'\n')
                                                #exit(1)
						if lineno >-1:
							line =line.replace(search, replace)
						sys.stderr.write(file+"\n")
                                                sys.stdout.write(line)

if len(sys.argv) < 3:
    print "usage: %s search_text replace_text directory extension" % os.path.basename(sys.argv[0])
    print "extension must have a leading dot, i.e. '.sel'"
    sys.exit(0)
search  = sys.argv[1]
replace = sys.argv[2]
path    = sys.argv[3]
exts    = sys.argv[4]
#print search, replace,path,exts 
searchreplace(path,search,replace,exts)
