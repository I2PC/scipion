# Makefile for Xmipp JNI libraries
CC = g++
SHAREDLIB = XmippJavaInterface.so
MAIN_OPTS = -s -pthread -shared
LIBS = -L/home/juanjo/xmipp/lib -lXmippData

SRC = xmipp_ImageDouble.cpp xmipp_MetaData.cpp xmipp_Projection.cpp
OBJS = $(SRC:.cpp=.os)
OPTS = -c -O2 -pthread -w -fPIC
INCLUDES = -I. -I/usr/lib/jvm/java-6-sun/include -I/usr/lib/jvm/java-6-sun/include/linux -I/home/juanjo/xmipp -I/home/juanjo/xmipp/libraries

JARFILE = XmippJavaInterface.jar
JAROPTS = -Mcf
PACKAGE = xmipp
JAVA_SRC = $(PACKAGE)/ImageDouble.java $(PACKAGE)/MDLabel.java $(PACKAGE)/MetaData.java $(PACKAGE)/Projection.java
JAVA_CLASSES = $(JAVA_SRC:.java=.class)
JAVAH_CLASSES = $(subst /,.,$(JAVA_SRC:.java=))

library_wrapper: java jar headers $(OBJS)
	g++ $(MAIN_OPTS) $(LIBS) $(OBJS) -o $(SHAREDLIB)

java:
	@echo " *** Compiling java sources..."
	./extract_labels.py
	javac $(JAVA_SRC)

jar:
	@echo " *** Generating jar file... $(JARFILE)"
	jar $(JAROPTS) $(JARFILE) $(JAVA_CLASSES)

headers:
	@echo " *** Generating header files..."
	javah $(JAVAH_CLASSES)

%.os:	%.cpp
	@echo " *** Compiling library: $@"
	$(CC) $? $(OPTS) $(INCLUDES) -o $@

clean:
	rm $(SHAREDLIB) $(OBJS) $(JAVA_CLASSES) $(JARFILE)

