# Spider Python Library
# Copyright (C) 2006, 2010 Health Research Inc.
#
# HEALTH RESEARCH INCORPORATED (HRI),
# ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455
#
# Email:  spider@wadsworth.org
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

import os, re, time
import struct, sys
from commands import getoutput
from types import *

WRITE_SPIREOUT = "WRITE_SPIREOUT"   # environmental variable for Spire

def fileReadLines(filename):
    "read a text file, return a list of lines"
    try:
        fp = open(filename,'r')
        B = fp.readlines()
        fp.close()
        return B
    except IOError, e:
        print 'Unable to open file \n' + filename, e
        return None

def fileWriteLines(filename, lines, append=0, mode='w'):
    "write a list of lines to a text file"
    mode = 'w'
    if append != 0 or mode != 'w':
        mode = 'a'
    fp = open(filename, mode)
    if type(lines) == type("string"):
        fp.write(lines)
    elif type(lines) == type(["list"]):
        fp.writelines(lines)
    fp.close()

###################################################################
#
# Reading and writing SPIDER document files
#
#     readdoc 
#     writedoc
#     nowisthetime, makeDocfileHeader, fixHeaders : functiones used by writedoc

# Usage: readdoc(file, column=2) or readdoc(file, columns=[2,3])
#
# Of they following 3 keywords (column, columns, lines), only one may be
# used in a call to readdoc. If none are used, then the default is "column=1"
#
# If column=int, returns that column as a list.
# If columns=[i,j,k], returns a tuple of lists.
# If lines=[i,j,k], returns a dictionary, each entry is a list of line values.
#
# Version 1.2 'keys' keyword is deprecated, use lines='all' instead
# If keys=0, returns a list of column values
#           (column 1 = 1st Spider doc file data column).
# If keys=1, returns a dictionary indexed by keys, uses column keyword.
# If keys='all', returns dictionary of all values, ignores column keyword.
# Only works for Spider Version 11.0 format, with spaces between columns.
def readdoc(filename, column=1, columns=None, line=None, lines=None, keys=0):
    "Read a SPIDER document file; return a list or a dictionary"
    # first, figure out the which keyword (columns, line, lines) is used
    columnlist = None
    linelist = None
    if lines != None:
        if isListorTuple(lines):
            linelist = lines
        elif lines == 'all':
            linelist = 'all'
        else:
            ci = checkInteger(lines)
            if ci != None:
                linelist = [ci]
            else:
                print "%s is not a valid value for 'lines' keyword" % str(lines)
                return None   # lines != None but can't parse it
    elif line != None:
        # try the 'line' keyword
        ci = checkInteger(line)
        if ci != None:
            linelist = [ci]
        elif isListorTuple(line):
            linelist = line
        else:
            print "%s is not a valid value for 'line' keyword" % str(line)
            return None
    elif columns != None:
        if isListorTuple(columns):
            columnlist = columns
        else:
            ci = checkInteger(columns)
            if ci != None:
                columnlist = [ci]
            else:
                return None   # columns != None but can't parse it
    else:
        # try the column keyword
        ci = checkInteger(column)
        if ci != None:
            columnlist = [ci]
        elif isListorTuple(column):
            columnlist = column
        else:
            return None

    # read the file data
    B = fileReadLines(filename)
    if B == None:
        return None

    # if keys=0, return a list of column values (or tuple of lists)
    if keys == 0 and columnlist != None:
        if len(columnlist) == 1:
            # only asking for 1 column
            if columnlist[0] != 0: col = columnlist[0] + 1      # first Spider data column
            else: col = 0      # key column
            C = []
            for line in B:
                d = line.split()
                try: int(d[0])
                except: continue
                
                try:
                    C.append(d[col])
                except:
                    ncols = int(d[1])
                    if col > ncols:
                        if len(C) > 0:
                            print "readdoc: column %d has no data at line %s" % (column,d[0])
                            del(B)
                            return map(float,C)
                        else:
                            print "column=%d, %s has only %d columns" % (column, filename, ncols)
                            del(B)
                            return None
            del(B)
            return map(float,C)
        else:
            # need to get multiple columns
            ncols = len(columnlist)
            listlist = []
            for c in columnlist:
                listlist.append([])  # a list of empty lists
            for line in B:
                d = line.split()
                try:  int(d[0])
                except: continue
                
                for i in range(ncols):
                    col = columnlist[i]
                    if col != 0: col = col + 1  # get past SPIDER's second column
                    try:
                        f = float(d[col])
                        listlist[i].append(f)
                    except:
                        pass
            del(B)
            return tuple(listlist)  # convert it to a tuple
        
    # return a dictionary of line values
    elif keys == 0 and linelist != None:
        D = {}
        if linelist == 'all':
            linedata = []
            for bline in B:
                s = bline.split()
                try: bkey = int(s[0])
                except: continue
                linedata =  map(float, s[2:])
                D[bkey] = linedata
            return D
        
        else: # a subset of lines
            for bline in B:
                s = bline.split()
                try: bkey = int(s[0])
                except: continue
                if bkey in linelist:
                    linedata =  map(float, s[2:])
                    D[bkey] = linedata

                    idx = linelist.index(bkey) # remove key from linelist
                    del(linelist[idx])
                    if len(linelist) == 0:
                        return D

    # if keys != 0, return a dictionary
    else:
        D = {}
        if keys == 'all':
            # each dictionary key points to list of all values in that line
            for line in B:
                d = line.split()
                try:  int(d[0])
                except: continue
                key = int(d[0])
                ncols = int(d[1]) + 1
                data = map(float,d[2:ncols+1])
                D[key] = data
        else:
            # dictionary contains a list of requested columns
            clist = []
            for col in columnlist:
                if col == 0: clist.append(col)
                else: clist.append(col+1)
            ncols = len(clist)
            for line in B:
                d = line.split()
                try:  int(d[0])
                except: continue
                key = int(d[0])
                data = []
                for col in clist:
                    data.append(float(d[col]))
                D[key] = data
        del(B)
        return D
# end readdoc -------------------------------------------------

# included for compatibility		
def readSpiderDocFile(filename, col_list=None):
    "Read a SPIDER document file; return a dictionary"
    return readdoc(filename, keys='all')

def numberOfColumns(docfile):
    " find the number of columns in a doc file"
    B = fileReadLines(docfile)
    for line in B:
     line = line.strip()
     if line == "" or line[0] == ";": continue

     d = line.split()
     if len(d) > 1:
         try:
             ncol = int(d[1])
             return ncol
         except:
             pass
    return 0

onespace = re.compile("\S+ \S+")

def combineHeaderSpaces(g,s,ncolumns):
    " g is list, s is string. Convert ['HAS', 'SPACE'] -> ['HAS SPACE']"
    m = onespace.search(s)
    # if an item with a space is found AND there are more headers than columns..
    if m and len(g) > ncolumns:
        start,end = m.span()
        substring = s[start:end]
        pp = substring.split()
        s1 = pp[0]
        s2 = pp[1]
        n = len(g)
        idx = -1
        for i in range(n-1):            
            if g[i] == s1 and g[i+1] == s2:
                idx = i
                break
        if idx > -1:
            # collapse list items i,i+1 to one item with a space
            h = g[:idx] + [g[idx] + ' ' + g[idx+1]]
            if n > idx+1:
                h = h + g[idx+2:]

            g = combineHeaderSpaces(h,s[end:],ncolumns)
            return g
        else:
            return g
    else:
        return g

def getDocfileHeaders(filename, output='list'):
    " tries to get the column headers from a doc file"
    try:
        fp = open(filename,'r')
    except:
        print 'unable to open %s' % (filename)
        return ""
    headers = []
    if output == "string":
        headers = ""
        
    max = 5
    i = 1
    while 1:
        s = fp.readline()
        if not s:
            break
        if len(s) > 0 and s[1] == ';':
            x = s.split()
            if len(x) == 5 and len(x[0]) == 8 and x[0][4] == '/' and x[2] == 'AT':
                # it's ;spi/dat   22-AUG-2008 AT 15:17:59   file.dat
                continue
            if len(x) > 1:
                if output == "string":
                    return s[2:].rstrip()
                
                g = x[1:]   # g is list, s is string

                # some headers have a space, so 2 items need to be combined
                ncol = numberOfColumns(filename)
                headers = combineHeaderSpaces(g,s.rstrip(), ncol)
                break
        i += 1
        if i > max:
            break

    fp.close()
    return headers


#########################################################################
# functions for writedoc

# returns e.g., ('16-OCT-03', '13:08:16', '031016130816')
def nowisthetime():
    "return current time as tuple of 3 strings: (date, time, ID)"
    tt = time.localtime(time.time())
    # localtime return format: (2003, 10, 16, 12, 48, 30, 3, 289, 1)
    #t = string.split(time.asctime(tt))
    t = time.asctime(tt).split()
    # asctime return format: 'Thu Oct 16 12:50:17 2003'
    mo = t[1].upper()
    day = t[2]
    if len(day) < 2: day = '0' + day
    timestr = t[3]
    yr = t[4]
    datestr = "%s-%s-%s" % (day, mo, yr)

    yr = yr[-2:]
    # this is just to get the month as a number
    d = map(str,tt)   # stringify all numbers in the tuple
    mon = d[1]
    if len(mon) < 2: mon = '0' + mon
    #(h,m,s) = string.split(timestr,':')
    (h,m,s) = timestr.split(':')
    idstr = "%s%s%s%s%s%s" % (yr,mon,day,h,m,s)

    return (datestr, timestr, idstr)

def makeDocfileHeader(filename, batext=None):
    "create the comment line used at the top of SPIDER document files"
    filename = os.path.basename(filename)
    fn, ext = os.path.splitext(filename)
    ext = ext[1:]
    if batext == None:
        batext = 'spl'   # Spider Python Library
    date,time,idstr = nowisthetime()
    h = " ;%s/%s   %s AT %s   %s\n" % (batext,ext,date,time,filename)
    return h

def fixHeaders(headers):
    "make all headers 11 characters in width; return doc string"
    w = 11
    docstr = " ; /    "
    for h in headers:
        d = len(h)
        if d > w:
            h = h[:w]
        docstr += h.rjust(w+1)
    docstr += "\n"
    return docstr

###################################################################
#
# type-checking functions

def checkInteger(x):
    " returns int if input can be converted to an integer, else None"
    try:
        i = int(x)  # works for 5, 5.0, '5'
        return i
    except:
        try:
            i = int(float(x))  # works for '5.0'
            return i
        except:
            return None
    return None

def isListorTuple(x):
    "returns 1 if input is a list or a tuple"
    if isinstance(x, ListType) or isinstance(x, TupleType) : return 1
    else: return 0
    
def isDictionary(d):
    "returns 1 if input is a Python dictionary"
    if isinstance(d, DictType): return 1
    else: return 0

def isListofLists(d):
    "returns 1 if input is a list, and 1st item in input is also a list"
    "actually works for tuples as well. Only checks 1st element "
    if not isListorTuple(d):
        return 0
    if len(d) < 1:
        return 0
    if isListorTuple(d[0]):
        return 1
    else:
        return 0
    
def getLastDocfileKey(docfile):
    "return the last key of a doc file"
    if not os.path.exists(docfile):
        return None
    cmd = 'tail %s' % docfile
    res = getoutput(cmd)
    s = res.split("\n")
    s.reverse()

    for line in s:
        if len(line) > 1 and line[1] != ";":
            ss = line.split()
            try:
                i = int(ss[0])
                return i
            except:
                pass
    return None
# --------------------------------------------------------------
# writedoc 
#    Data can be organized as columns or line. Call to writedoc should
#    use EITHER columns OR lines.
#    Data must be integer or float. (they can be in string format)
#    columns: a list of lists; each doc file column is a list
#    lines : a list of lists; each doc file line is a list (w/o key)
#    headers: a list of strings
#
# todo: check if columns have different lengths
def writedoc(filename, columns=None, lines=None, headers=None, keys=None, mode='w'):
    "write data to a file in SPIDER document file format"
    if not isListofLists(columns) and not isListofLists(lines):
        if isDictionary(columns):
            return writeSpiderDocFile(filename, columns, headers=headers, mode=mode)
        else:
            print "writedoc: columns or lines must be a list of lists"
            return

    "filename must have data extension"
    try:
        fp = open(filename, mode)
    except:
        print "Unable to open %s for writing." % filename
        return

    # write Spider doc file header
    lastkey = None
    _APPEND = 0
    if mode == 'w':
        hdr = makeDocfileHeader(os.path.basename(filename))
        fp.write(hdr)
    elif mode == 'a':
        _APPEND = 1
        try:
            lastkey = getLastDocfileKey(filename)
        except:
            pass
    # write column headings
    if headers != None and type(headers) == type(["list"]):
        fp.write(fixHeaders(headers))

    datalines = []

    # write data columns
    if columns != None:
        ncol = len(columns)  # number of columns
        n = len(columns[0])  # length of 1st column (assumes all have same length)
        if keys == None:
            if lastkey == None:
                keys = range(1,n+1)
            else:
                keys = range(lastkey+1, lastkey+n+1)

        for i in range(n):
            dstr = "%5d %2d" % (int(keys[i]), int(ncol))
            for j in range(ncol):
                dstr += " %11g" % float(columns[j][i])
            datalines.append(dstr+"\n")

    # write data lines
    elif lines != None:
        n = len(lines)       # number of lines
        if keys == None:
            if lastkey == None:
                keys = range(1,n+1)
            else:
                keys = range(lastkey+1, lastkey+n+1)

        for i in range(n):
            line = lines[i]
            ncol = len(line) # number of columns
            dstr = "%5d %2d" % (int(keys[i]), ncol)
            for item in line:
                dstr += " %11g" % float(item)
            datalines.append(dstr+"\n")

    fp.writelines(datalines)
    fp.close()

    if not _APPEND: # i.e. it's a new doc file
        writeSpireoutFile(filename)

#included for compatibiity:
# data must be in dictionary form: D[key] = [list of column values]
def writeSpiderDocFile(filename, data, headers=None, append=0, mode='w'):
    "write data (in dictionary form) to a file in SPIDER document file format"
    if append > 0: mode = 'a'
    if mode == 'a': _APPEND = 1
    else: _APPEND = 0
    
    if not isDictionary(data):
        # if it's not a dictionary, see if it's a list of lists
        if isListofLists(data):
            return writedoc(filename, columns=data, headers=headers, mode=mode)
        else:
            return 0
    try:
        fp = open(filename, mode)
    except:
      print "unable to open %s for writing" % filename
      return 0
    # write Spider doc file header
    hdr = makeDocfileHeader(os.path.basename(filename))
    fp.write(hdr)
    # and any column headings
    if headers != None and type(headers) == type(["list"]):
        fp.write(fixHeaders(headers))

    # write data
    keys = data.keys()
    keys.sort()
    if len(keys) > 0:
        firstkey = keys[0]
        if firstkey in data:
            v1 = data[firstkey]
        else:
            print "writeSpiderDocFile: key not found in data"
    else:
        print "writeSpiderDocFile: empty data dictionary"
        return 0
        
    if isinstance(v1, ListType):
        for key in keys:
            values = data[key]
            n = len(values)
            h = "%5d %2d " % (int(key),int(n))
            for value in values:
                h += " %11g " % (float(value))
            fp.write(h+"\n")
    else:
        # it's supposed to be a list! But if it's not..
        for key in keys:
            value = data[key]
            try:
                f = float(value)
            except:
                print "writedoc: unable to convert %s" % str(value)
                print "writedoc: Dictionary elements must be lists!"
                return 0
            h = "%5d %2d %11g\n" % (int(key), 1, f)
            fp.write(h)
        
    fp.close()
    if not _APPEND: # i.e. it's a new doc file
        writeSpireoutFile(filename)
    return 1

def writeSpireoutFile(filename):
    # this section added Aug 2006, for Spire compatibility
    if WRITE_SPIREOUT in os.environ:
        spireout = os.environ[WRITE_SPIREOUT]
        if os.path.exists(spireout):
            fp = open(spireout, mode='a')
        else:
            fp = open(spireout, mode='w')
        date, time, id = nowisthetime()
        fn = os.path.basename(filename)
        s = "  %s AT %s    OPENED NEW DOC FILE: %s\n" % (date, time, fn)
        fp.write(s)
        fp.close()

###################################################################
# convenience functions
#
def list2int(a):
    " converts a list of strings to integers"
    return map(int, map(float,a))

def list2float(a):
    " converts a list of strings to floats"
    return map(float,a)

# given ("mic****", 89) returns "mic0089"
# Substitutes the first set of asterisks it finds (i.e. leftmost).
# If number of *'s too small for number, filename is extended.
# asterisks can be replaced by another char
def makeFilename(filename, n, char='*'):
    "substitutes asterisks for number: ('mic****', 89) returns 'mic0089'"
    try:
        n = int(float(n))
    except:
        print "makeFilename: unable to convert %s to integer" % str(n)
    number = str(n)

    re_findchar = re.compile('[%s]+' % char)

    a = re_findchar.search(filename)
    if not a:
        print "makeFilename: no '%s' found in %s" % (char, filename)
        return ""

    (start,end) =  a.span()
    sp = (end - start)
    num = number.zfill(sp)  # pad the numeric string with zeroes out to req. length
    f = filename[:start] + num + filename[end:]
    return f

re_asterisk = re.compile('[*]+')

# old version of makeFilename included for compatibility
def makeSpiderFilename(filename, n):
    "substitutes asterisks for number: ('mic****', 89) returns 'mic0089'"
    try:
        n = int(float(n))
    except:
        print "makeSpiderFilename: unable to convert %s to integer" % str(n)
    number = str(n)

    a = re_asterisk.search(filename)
    if not a:
        #print "makeSpiderFilename: no asterisks in %s" % filename
        return filename

    (start,end) =  a.span()
    sp = (end - start)
    num = number.zfill(sp)  # pad the numeric string with zeroes out to req. length
    f = filename[:start] + num + filename[end:]
    return f


###################################################################
# File number routines:
#   filenumber(filename) returns an integer
#   getfilenumber(filename) returns string with leading zeroes
#   name2template(filename) given 'mic021.dat', returns 'mic***.dat'
#   template2filename(template, n) given (pic***.dat, 3) returns pic003.dat
#   numberlist2string(nlist) converts [1,2,3,4,6,8] -> '1-4,6,8'
#   range2list(nrange) given string '1-4', outputs list [1,2,3,4]

#re_nums = re.compile('\d+\D?')  # ints followed by one non-int char

def filenumber(file):
    "returns file number (integer nearest the file extension)"
    if len(file) == 0: return None
    n = getfilenumber(file)
    if n:
        return int(n)
    else:
        return None

def getfilenumber(filename):
    "returns file number as a string with leading zeroes "
    filename = os.path.basename(filename)
    fname,ext = os.path.splitext(filename)

    numstr = ""
    f = list(fname)
    f.reverse()
    done = 0
    for ch in f:
        if not done:
            try:
                int(ch)
                numstr = ch + numstr
            except:
                if numstr != "":
                    done = 1
    return numstr
            
def name2template(filename, all=0):
    " given 'mic021.dat' --> returns mic***.dat "
    " by default, only replaces number nearest extension. all !=0 replaces all"
    if len(filename) == 0: return ""
    path, basename = os.path.split(filename)
    fname,ext = os.path.splitext(basename)
    
    newfn = ""
    f = list(fname)
    f.reverse()
    if all:
        for ch in f:
            try:
                int(ch)
                newch = '*'
            except:
                newch = ch
            newfn = newch + newfn
    else:
        found = 0
        for ch in f:
            if not found:
                try:
                    int(ch)
                    newch = '*'
                except:
                    newch = ch
                    if newfn and newfn[0] == '*':
                        found = 1
            else:
                newch = ch
            newfn = newch + newfn
    fname = os.path.join(path,newfn) + ext
    return fname

# template should have asterisks, num can be int or a numbered filename.
# Like makeSpiderFilename, but it can accept a filename instead of a number.
def template2filename(template, n=0):  #numfile=None, n=None):
    "replaces asterisks with number: (pic***.dat, doc003.dat) returns pic003.dat"
    if type(n) == type(1):
        pass
    elif type(n) == type("string"):
        n = filenumber(n)
    else:
        print "template2filename: unable to parse input"
        return ""
    nstars = template.count("*")
    if nstars == 0:
        return template
    if len(str(n)) > nstars:
        print "template2filename: **** Warning number larger than template"
    numstr = str(n).zfill(nstars)
    sts = "*" * nstars
    filename = template.replace(sts,numstr)
    return filename

# list2range: hyphenates runs of consecutive intgers.
#     input: list [1,2,3,6,7,8,9,10,11,17]
#     output: list with strings ['1-3', '6-11', '17']
def list2range(f):
    " input is a list of integers "
    fn = []
    for item in f:
        try:
            fn.append(int(item))
        except:
            print "list2range: unable to convert %s to integer" % str(item)
            return []
    if len(fn) < 1: return []
    fn.sort()
    N = []
    p = fn[0]
    start = p
    n = p
    fn = fn[1:]

    while len(fn) > 0:
        n = fn[0]
        if n == p+1:
            p = n
        elif n == p:
            print "list2range: Warning, %s occurs more than once" % str(n)
            p = n
        else:
            if start != p:
                if p == start+1: # don't create two-item ranges (eg, '1-2')
                    N.append(str(start)) ; N.append(str(p))
                else:
                    next = "%s-%s" % (str(start),str(p))
                    N.append(next)
            else:
                N.append(str(start))
            start = n
            p = n
        fn = fn[1:]

    if start != n:
        if n == start+1:
            N.append(str(start)) ; N.append(str(n))
        else:
            next = "%s-%s" % (str(start),str(n))
            N.append(next)
    else:
        N.append(str(start))

    return N

# input: list of strings (output from list2range)
# output: string of concatenated items
def range2string(n): 
    " converts ['1-4','7'] to '1-4,7' "
    if n == None or len(n) == 0:
        return ""
    s = ""
    for item in n:
        s += item + ','
    s = s[:-1] # remove trailing comma
    return s

def numberlist2string(numberlist):
    " a list of integers is converted to a string, "
    " hyphenating consecutive numbers. [1,2,3,4] -> '1-4' "
    d = list2range(numberlist)
    return range2string(d)

def range2list(numberstring):
    "input string: '1-4' -> output list: [1,2,3,4] "
    if numberstring == "" or None:
        return []
    L = numberstring.split(',')
    K = []
    for item in L:
        if item.find('-') > -1:
            xlist = item.split('-')
            start = int(xlist[0])
            end = int(xlist[-1])
            for i in range(start,end+1):
                K.append(i)
        else:
            K.append(int(item))
    return K

###################################################################
#
# Checking file types
#
#     istextfile()       boolean
#     isSpiderDocFile()  boolean
#     isSpiderImage()    boolean
#     isSpiderBin()      returns "image","volume","Fourier"  or 0

# ------ text file functions -------------

noNumbers = re.compile("[^\d^\s^\.^\+^\-Ee]")
text_characters = "".join(map(chr, range(32, 127)) + list("\n\r\t\b"))
from string import maketrans
_null_trans = maketrans("", "")

# Applications importing this file should call istextfile, not istext,
# which is used internally.
# istextfile returns 1 for text, 0 for binary, 0 for error (not found?)
# pdf test added, cos they can get either answer.
def istextfile(filename, blocksize = 512):
    "returns 1 if input is a text file (pdf's and zero-length files are binary)" 
    if os.path.isdir(filename):
        return 0
    name,ext = os.path.splitext(filename)
    if ext.lower() == ".pdf":
        return 0
    try:
        res = istext(open(filename).read(blocksize))
        return res
    except:
        return 0

def istext(s):
    "returns 1 if input is a text file (pdf's and zero-length files are binary)" 
    if "\0" in s:
        return 0

    if not s:  # Empty files 
        return 1

    # Get the non-text characters (maps a character to itself then
    # use the 'remove' option to get rid of the text characters.)
    t = s.translate(_null_trans, text_characters)

    # If more than 30% non-text characters, then
    # this is considered a binary file
    ratio = float(len(t)) / float(len(s))
    if ratio > 0.30:
        return 0
    return 1

# Quits as soon as it gets a good data line, i.e.,
# int1 int2 [floats], where no. floats = int2
def isSpiderDocfile(file):
    "returns 1 if input is a SPIDER document file"
    try:
        fp = open(file, 'r')
    except:
        print 'unable to open %s' % (file)
        return 0

    comments = 0
    isDoc = 0
    blank = 0
    while 1:
        s = fp.readline()
        if s == "":  # only EOF should return blank
            break

        if len(s) > 2 and s[0] == " " and s[1] == ';':   # Spider comment line
            continue

        if noNumbers.match(s):  # if find any nondigits, +, _ etc
            isDoc = 0
            break

        ss = s.split()
        # test for new format: nums divided by blanks, 1st value is an int,
        try:
            i = int(ss[0])
            # and there are N data columns, where N = s[1]
            n = int(ss[1])
            if len(ss[2:]) >= n:
                try:
                    float(ss[2])  # we'll just test one
                    isDoc = 1
                except:
                    isDoc = 0
                break         # then it's new (SPIDER 11.0 Feb 2004)
        except:
            pass
        
        # see if it's the older fixed column format
        if len(s) < 13:
            isDoc = 0
            break
        try:
            key = int(s[0:6])   # 1st 6 chars are key
            n = int(s[6])       # 7th char is N
            f = float(s[7:13])   # see if there's 1 good data value
            isDoc = 1
            break
        except:
            isDoc = 0
            break
    fp.close()
    return isDoc

def stripComment(line, strip=1):
    "removes all text after and including the 1st semicolon in a string"
    n = line.find(";")
    if n > -1:
        line = line[:n].rstrip()
        
    if strip:
        line = line.strip()
    return line

re_hdr = re.compile('END +BATCH +HEADER')  
re_reg = re.compile('[xX][0-9][0-9] *=')      # "x11 =" patterns
re_nam = re.compile('\[[ a-zA-Z0-9_-]+\] *=') # "[symbol] =" patterns
re_sym = re.compile('\?[ \w]+\?')   # ?text? patterns
# top line of procedures
re_reglist = re.compile('([xX][0-9][0-9])(, *[xX][0-9][0-9])*') # x11,x12,x13
re_namlist = re.compile('(\[[ a-zA-Z0-9_-]+\])(, *\[[ a-zA-Z0-9_-]+\])*')
re_nam1 = re.compile('(\[[ a-zA-Z0-9_-]+\])')  # named reg in proc hdr

# This is not really dependable - there are just too many variants to
# catch them all. Plus it may return with false positives.
def isSpiderBatchfile(file):
    "returns 1 if input is a SPIDER batch file"
    """ only checks first few lines of text.
        Returns 1 if finds any of the following:
            ; --- End batch header ---
            "FR G/L" pattern followed by [symbol] on next line
            "x11=" register assignment pattern
            "[symbol]=" named register assignment
        Returns 2 if it thinks its a procedure, with the first line:
            [x11,x12]
            ([ang-step],[ang-limit],[radius])
            again, there are too many variants to catch them all
    """
    #B = fileReadLines(file)        takes too long for huge text files
    #if B == None or len(B) == 0:
    #   return 0

    max = 40
    B =[]
    fp = open(file,'r')
    for i in range(max):
        try:
            B.append(fp.readline())
        except:
            pass
    fp.close()

    nlines = len(B)
    if nlines < 40:
        max = nlines

    for i in range(max):
        line = B[i]
        if not line:
            return 0
        if line.find("RESULTS FILE FLUSHED") > -1:
            return 0
        line = line.strip()
        line = line.upper()
        if len(line) < 2: continue

        # test if 1st line is a procedure call
        if (line[0]=='[' and line[-1]=="]") or (line[0]=='(' and line[-1]==")"):
            if re_reglist.match(line[1:-1]):
                return 2
            if re_namlist.match(line[1:-1]):
                return 2
        
        cmd = ""
        if len(line) > 3:
            cmd = line[0:4]
            
        if len(line) == 0:
            continue
        elif re_hdr.search(line):
            #print "hdr: " + line
            return 1
        # comment check must come after header check
        elif line[0] == ";":
            if line.find('SPIDER') > -1:  # a comment with the word 'spider'?
                return 1
            else:
                continue
        elif re_reg.match(line):
            #print "reg: " + line
            return 1
        elif re_nam.match(line):
            #print "nam: " + line
            return 1
        elif re_sym.match(line):
            #print "sym: " + line
            return 1
        elif cmd == "FR G" or cmd == "FR L":
            #print "FR: " + line
            return 1
        elif line == "FR":
            nextline = B[i+1].strip()
            if re_sym.match(nextline):
                return 2

    return 0

def isSpiderProcedurefile(file):
    "returns 1 if input is a SPIDER procedure file"
    if isSpiderBatchfile(file) == 2:
        return 1
    else:
        return 0
    
# ------ binary file functions -------------

def isInt(f):
    "returns 1 if input is an integer"
    try:
        i = int(f)
        if f-i == 0: return 1
        else:        return 0
    except:
        return 0

iforms = [1,3,-11,-12,-21,-22]

# returns header tuple, if t is a valid Spider header,
# otherwise returns 0
def isSpiderHeader(t):
    "returns tuple of values from a valid SPIDER header, else 0"
    h = (99,) + t   # add 1 value so can use spider header index start=1
    # header values 1,2,5,12,13,22,23 should be integers
    for i in [1,2,5,12,13,22,23]:
        if not isInt(h[i]): return 0
    # check iform
    iform = int(h[5])
    if not iform in iforms: return 0
    # check other header values
    labrec = int(h[13])   # no. records in file header
    labbyt = int(h[22])   # total no. of bytes in header
    lenbyt = int(h[23])   # record length in bytes
    #print "labrec = %d, labbyt = %d, lenbyt = %d" % (labrec,labbyt,lenbyt)
    if labbyt != (labrec * lenbyt): return 0
    # looks like a valid header
    return h

# returns "image","volume","Fourier"  or 0
def isSpiderBin(filename):
    "returns nonzero value if input is a SPIDER binary file"
    if not os.path.exists(filename):
        return 0
    minsize =  27 * 4  # 27 floating points
    if os.path.getsize(filename) < minsize:
        return 0
    try:
        fp = open(filename,'rb')
        f = fp.read(minsize)   # read 27 * 4 bytes
        fp.close()
    except:
        return 0
    bigendian = 1
    t = struct.unpack('>27f',f)    # try big-endian first
    hdr = isSpiderHeader(t)
    if hdr == 0:
        bigendian = 0
        t = struct.unpack('<27f',f)  # little-endian
        hdr = isSpiderHeader(t)
    if hdr == 0:
        return 0
    
    iform = int(hdr[5])
    if iform == 1:
        istack = hdr[24]
        if istack == 0:
            return "image"
        else:
            return "stack"
    elif iform == 3:
        return "volume"
    elif iform in [-11,-12,-21,-22]:
        return "Fourier"
    else:
        return 0

def isSpiderImage(file):
    "returns 1 if input is a SPIDER 2D image"
    if isSpiderBin(file) == "image": return 1
    else: return 0

def isSpiderVolume(file):
    "returns 1 if input is a SPIDER 3D volume"
    if isSpiderBin(file) == "volume": return 1
    else: return 0

def isSpiderStack(file):
    "returns 1 if input is a SPIDER stack file"
    if isSpiderBin(file) == "stack": return 1
    else: return 0
    
###################################################################
#
# Utilities for finding and testing SPIDER

def testSpider(spider):
    "returns 1 if input is a working path to SPIDER"
    file = 'test6637'
    ext = ".bat"
    filename = file + ext
    fp = open(filename, 'w')
    fp.write("en d\n")
    fp.close()
    spicmd = "%s bat/dat @%s" % (spider, file)

    success = 0
    output = getoutput(spicmd)
    if output.find('Results file') > 0:
        success = 1
        log = "LOG" + ext
        if os.path.exists(log):
            os.remove(log)
    os.remove(filename)
    return success

def programExists(prog):
    "a wrapper for os.path.exists that won't crash"
    try:
        if os.path.exists(prog): return 1
        else: return 0
    except:
        return 0

def findProgram(prog):
    "Use the Unix 'which' command to find a program"
    if os.name != 'posix':
        print 'not a posix system: no "which" command?'

    cmd = 'which %s' % prog
    out = getoutput(cmd)
    # output from 'which' command may contain newlines and spaces
    if out.find(os.linesep) > -1:
        lines = out.split(os.linesep)
        for line in lines:
            if line.find(" ") > -1:
                d = line.split()
                for item in d:
                    if programExists(item):
                        return item
            else:
                if programExists(line):
                    return line
    elif out.find(" ") > -1:
        d = out.split()
        print d
        for item in d:
            if programExists(item):
                return item
    else:
        if programExists(out):
            return out
    # failure
    return ""

def findSpider():
    "returns path to SPIDER, or else an empty string"
    spider = findProgram('spider')
    if spider != "" and testSpider(spider):
        return spider
    else:
        return ""

def runSpider(spider, batch, dataext):
    bat, batext = os.path.splitext(batch)
    batext = batext[1:]
    if dataext[0] == '.':
        dataext = dataext[1:]
    cmd = "%s %s/%s @%s" % (spider, batext, dataext, bat)
    out = getoutput(cmd)
    return out

###################################################################
#
# Reading and writing the SPIDER header
#
SpiderHeaderDict = { 1 : 'nslice ', 2 : 'nrow ', 3 : 'irec ', 4 : 'nhistrec ',
                     5 : 'iform ', 6 : 'imami ', 7 : 'fmax ', 8 : 'fmin ',
                     9 : 'av ', 10 : 'sig ', 11 : 'ihist ', 12 : 'nsam ',
                     13 : 'labrec ', 14 : 'iangle ', 15 : 'phi ', 16 : 'theta ',
                     17 : 'gamma ',18 : 'xoff ',19 : 'yoff ',20 : 'zoff ',
                     21 : 'scale ', 22 : 'labbyt ', 23 : 'lenbyt ', 24 : 'istack ',
                     25 : 'NOTUSED ', 26 : 'maxim ', 27 : 'imgnum ', 28 : 'lastindx ',
                     29 : 'unused ', 30 : 'unused ', 31 : 'Kangle ', 32 : 'phi1 ',
                     33 : 'theta1 ', 34 : 'psi1 ', 35 : 'phi2 ', 36 : 'theta2 ',
                     37 : 'psi2 '}

# item[0] (bigendian flag) is not part of the Spider header. It is added
# so that SPIDER indices (starting with 1) may be used.
# hdr is the array returned by getSpiderHeader.
class SpiderHeaderClass:
    def __init__(self, hdr):
        self.header = hdr
        self.hdrlen = len(hdr)
        self.bigendian = hdr[0]

        for i in range(1, self.hdrlen):
            if i in SpiderHeaderDict:
                name = SpiderHeaderDict[i]
                if name in ['NOTUSED', 'unused']:
                   continue
                val = hdr[i]
                s = "self.%s = %f" % (name, val)
                exec(s)
        
        if self.hdrlen > 9:
            self.avg = hdr[9]  # alternate access format
        if self.hdrlen > 31:
            self.kangle  = hdr[31]

# Create a SPIDER header for binary files
def makeSpiderHeader(dims):
    " dims must be (nsam, nrow), or (nsam, nrow, nslice) "
    if len(dims) == 2:
        nsam, nrow = dims[0], dims[1]
        nslice = 1.0
        iform = 1.0
        isVolume = 0
    elif len(dims) == 3:
        nsam, nrow, nslice = dims[0], dims[1], dims[2]
        iform = 3.0
        isVolume = 1
    else:
        return []

    lenbyt = nsam * 4  # There are labrec records in the header
    labrec = 1024 / lenbyt
    if 1024%lenbyt != 0: labrec += 1
    labbyt = labrec * lenbyt
    hdr = []
    nvalues = labbyt / 4
    for i in range(nvalues):
        hdr.append(0.0)
        
    if len(hdr) < 23:
        return []

    # NB these are Fortran indices
    hdr[1]  = float(nslice) # nslice (=1 for an image) 
    hdr[2]  = float(nrow)   # number of rows per slice
    hdr[5]  = iform         # iform for 2D image
    hdr[12] = float(nsam)   # number of pixels per line
    hdr[13] = float(labrec) # number of records in file header
    hdr[22] = float(labbyt) # total number of bytes in header
    hdr[23] = float(lenbyt) # record length in bytes

    # adjust for Fortran indexing
    hdr = hdr[1:]
    hdr.append(0.0)
    # pack binary data into a string
    hdrstr = []
    for v in hdr:
        hdrstr.append(struct.pack('f',v))
    return hdrstr

def getSpiderHeader(filename, n=27):
    " returns first n numbers, with Spider indices (starting at 1)"
    " if n = 'all', returns entire header "
    if not os.path.exists(filename):
        return 0
    getall = 0
    if not isInt(n):
        n = 27
        getall = 1
    nwords = n * 4  # no. floating point words 
       
    if os.path.getsize(filename) < nwords:
        return 0
    try:
        fp = open(filename,'rb')
        f = fp.read(nwords)   # read 27 * 4 bytes
        fp.close()
    except:
        return 0
    bigendian = 1
    bigformat = '>%df' % n
    t = struct.unpack(bigformat,f)    # try big-endian first
    hdr = isSpiderHeader(t)
    if hdr == 0:
        bigendian = 0
        littleformat = "<%df" % n
        t = struct.unpack(littleformat,f)  # little-endian
        hdr = isSpiderHeader(t)

    if hdr == 0:
        return 0
    else:
        # check if user requested the entire header
        if getall:
            labbyt = int(hdr[22])   # total no. of bytes in header
            hdr = getSpiderHeader(filename, n=labbyt)
        hdr = list(hdr)
        hdr[0] = bigendian
        return hdr

    
# returns [type,  (dimensions), and if there are any,(stats) ]
# where type = "image","volume","Fourier","stack" (but only for image stacks)
# or returns 0
def spiderInfo(filename):
    if not os.path.exists(filename):
        return 0

    type = isSpiderBin(filename)
    if type == 0:
        return 0
    hdr = getSpiderHeader(filename)  # header with Spider indices (starting at 1)
    info = getSpiderInfo(hdr)
    return [type] + info

# return [ (dimensions), (stats) ]  
def getSpiderInfo(h):
    " assumes its a valid header "
    #h = (99,) + t   
    nsam = int(h[12])
    nrow = int(h[2])
    nslice = int(h[1])
    iform = int(h[5])

    dim2D = [1, -11, -12]
    dim3D = [3, -21, -22]
    if iform in dim3D:
        dims = (nsam, nrow, nslice)
    else:
        dims = (nsam, nrow)
        
    imami = int(h[6])
    if imami != 0:
        max = float(h[7])
        min = float(h[8])
        avg = float(h[9])
        std = float(h[10])
        stats = (max, min, avg, std)
        return [dims, stats]
    else:
        return [dims]
   
