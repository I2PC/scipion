#!/usr/bin/env python
FILENAMENUMBERLENTGH=6
#---------------------------------------------------------------------------
# utils_xmipp.composeFileName (Xmipp-like)
#---------------------------------------------------------------------------
def composeFileName(rootname,number,extension):
    if (number != -1):
        output = rootname + str(number).zfill(FILENAMENUMBERLENTGH)
    else:
        output = rootname

    if (extension != ''):
        output += '.' + extension 

    return output

def composeWildcardFileName(rootname,extension):
    output = rootname
    for i in range(FILENAMENUMBERLENTGH):
        output += '?'

    if (extension != ''):
        output += '.' + extension 

    return output

def getCommaSeparatedIntegerList(inputstring):
    import string
    lista=string.split(inputstring,",")
    output=[]
    for i in range(len(lista)):
        interval=string.split(lista[i],'-')
        if len(interval)==1:
            if not interval[0]=='':
                output+=[int(interval[0])]
        else:
            output+=range(int(interval[0]),
                        int(interval[1])+1)
    return output

import os
def unique_filename(file_name):
    counter = 1
    file_name_parts = os.path.splitext(file_name) # returns ('/path/file', '.ext')
    while os.path.isfile(file_name): 
        file_name = file_name_parts[0] + '_' + str(counter) + file_name_parts[1]
        counter += 1
    return file_name

