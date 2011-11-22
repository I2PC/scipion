# Makefile for Xmipp JNI libraries
XMIPP_DIR = $$XMIPP_HOME

CC = g++
SHAREDLIB = XmippJavaInterface.so
MAIN_OPTS = -s -pthread -shared
LIBS = -L$(XMIPP_DIR)/lib -lXmippData

SRC = xmipp_ImageDouble.cpp xmipp_ImageGeneric.cpp xmipp_ImageGeneric_.cpp xmipp_MetaData.cpp xmipp_Projection.cpp xmipp_ExceptionsHandler.cpp xmipp_Aux.cpp
OBJS = $(SRC:.cpp=.os)
OPTS = -c -O2 -pthread -w -fPIC
INCLUDES = -I. -I/usr/lib/jvm/java-6-sun/include -I/usr/lib/jvm/java-6-sun/include/linux -I$(XMIPP_DIR) -I$(XMIPP_DIR)/libraries

JARFILE = XmippJavaInterface.jar
JAROPTS = -Mcf
SRCDIR = src/
PACKAGE = xmipp
IJ_JAR = $(XMIPP_DIR)/external/imagej/ij.jar
JAVA_FILES = ImageDouble.java ImageGeneric.java ImageGeneric_.java Filename.java \
MDLabel.java MetaData.java Projection.java CTFDescription.java Program.java
JAVA_SRC = $(JAVA_FILES:.java=.class)

$(SRCDIR)$(PACKAGE)/ImageDouble.java $(SRCDIR)$(PACKAGE)/ImageGeneric.java \
$(SRCDIR)$(PACKAGE)/ImageGeneric_.java $(SRCDIR)$(PACKAGE)/Filename.java \
$(SRCDIR)$(PACKAGE)/MDLabel.java $(SRCDIR)$(PACKAGE)/MetaData.java \
$(SRCDIR)$(PACKAGE)/Projection.java $(SRCDIR)$(PACKAGE)/CTFDescription.java \
$(SRCDIR)$(PACKAGE)/Program.java
JAVA_CLASSES = $(JAVA_SRC:.java=.class)
JAVAH_CLASSES = $(subst /,.,$(JAVA_SRC:.java=))

TESTS_PACKAGE = tests
TESTS_SRC = $(TESTS_PACKAGE)/ImageTest.java $(TESTS_PACKAGE)/MetaDataTest.java $(TESTS_PACKAGE)/ProjectionTest.java $(TESTS_PACKAGE)/PreviewTest.java
TESTS_CLASSES = $(TESTS_SRC:.java=.class)

library_wrapper: java jar headers $(OBJS) javatests
	@echo " *** Compiling main lib... $(SHAREDLIB)"
	g++ $(MAIN_OPTS) $(LIBS) $(OBJS) -o $(SHAREDLIB)

java:
	@echo " *** Compiling java sources..."
#	./extract_labels.py
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

javatests: $(TESTS_SRC)
	@echo " *** Compiling tests"
	javac -cp $(IJ_JAR):./* $(TESTS_SRC)

clean:
	rm $(SHAREDLIB) $(OBJS) $(JAVA_CLASSES) $(JARFILE) $(TESTS_CLASSES)

