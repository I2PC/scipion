#!/usr/bin/python

import os, sys;
    
f = open("../../data/metadata_label.h");    
fOut = open("xmipp/MDLabel.java", 'w+')

counter = 0;
fOut.write("package xmipp; \n")
fOut.write("public class MDLabel {\n")

for line in f:
    l = line.strip();
    if l.startswith("MDL_"):
        if l.startswith("MDL_LAST_LABEL"):
            l = l.replace("MDL_LAST_LABEL", "MDL_LAST_LABEL = " + str(counter) + ";")   
        if (l.find("=") == -1):
            l = l.replace(",", " = " + str(counter) + "; ")   
            counter = counter + 1;
        else:
            l = l.replace(",", ";")                                       
            
        fOut.write("   public static final int %s\n" % l)
        
fOut.write("}\n")
    
