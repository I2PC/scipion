#!/usr/bin/env python
import pickle
import numpy
import re

if __name__ == '__main__':
    vec = open( "vec_ani.txt", "r" )
    data= []
    line_number=0
    for aLine in vec:
        fields= re.split("[ ]+",aLine)
        line_number+=1
        #print line_number, fields[0], fields[1]
        data.append( fields )
    vec.close()
    f = open( "vec_ani.pkl", "w" )
    pickle.dump(data, f)
    f.close()




