CC=g++
CFLAGS=-Wall -fpic
LIBS=-I"/usr/lib/jvm/java-6-sun-1.6.0.20/include/" -I"/usr/lib/jvm/java-6-sun-1.6.0.20/include/linux/"
LDFLAGS=-shared
OUTPUTLIB=libXmippDataJava.so
SWIG_INTERFACE=tmp.i#swigXmippDataJava.i
WRAPPER=tmp_wrap.cxx#swigXmippDataJava_wrap.cxx

all: swig $(OUTPUTLIB) java

$(OUTPUTLIB):
	@echo "Compiling shared library $(OUTPUTLIB)"
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(OUTPUTLIB) $(LIBS)

swigXmippDataJava.o:	classes_wrap.cxx
	@echo "Compiling $@"
	$(CC) $(CFLAGS) -c classes_wrap.cxx $(LIBS)

clean:
	rm -f *.o *~ *.class classes_wrap.cxx classes.java Point2D.java Point3D.java classesJNI.java $(OUTPUTLIB)

swig:
	@echo "SWIG: generating java files"
	swig -java -c++ $(SWIG_INTERFACE)

java:
	javac *.java

