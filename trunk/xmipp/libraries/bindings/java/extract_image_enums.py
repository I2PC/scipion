#!/usr/bin/python

import os, sys;

def extract_enum(header_file,java_file,class_name,pattern):
	f = open(header_file);    
	fOut = open(java_file, 'w+')
	
	counter = 0;
	fOut.write("package xmipp; \n")
	fOut.write("public class " + class_name + " {\n")
	for line in f:
	    l = line.strip();
	    if l.startswith(pattern):
	    	last_label_pattern = pattern + "LAST_LABEL"
	        if l.startswith(last_label_pattern):
	            l = l.replace(last_label_pattern, last_label_pattern +" = " + str(counter) + ";")   
	        if (l.find("=") == -1):
	            l = l.replace(",", " = " + str(counter) + "; ")   
	            counter = counter + 1;
	        else:
	            l = l.replace(",", ";")                                       
	            
	        fOut.write("   public static final int %s\n" % l)    

	fOut.write("}\n")
    
extract_enum("../../data/image_base.h","xmipp/ImageWriteMode.java","ImageWriteMode","WRITE_")
extract_enum("../../data/image_base.h","xmipp/CastWriteMode.java","CastWriteMode","CW_")