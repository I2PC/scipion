echo "Compiling java sources..."
./extract_labels.py
javac xmipp/*.java
jar -Mcf XmippJavaInterface.jar xmipp

echo "Generating header files"
javah xmipp.ImageDouble xmipp.Projection xmipp.MetaData

echo "Compiling library"
g++ -o ImageDouble.os -c -O2 -pthread -w -O2 -pthread -fPIC -I. -I/usr/lib/jvm/java-6-sun/include -I/usr/lib/jvm/java-6-sun/include/linux  -I/home/juanjo/xmipp -I/home/juanjo/xmipp/libraries xmipp_ImageDouble.cpp

g++ -o MetaData.os -c -O2 -pthread -w -O2 -pthread -fPIC -I. -I/usr/lib/jvm/java-6-sun/include -I/usr/lib/jvm/java-6-sun/include/linux  -I/home/juanjo/xmipp -I/home/juanjo/xmipp/libraries xmipp_MetaData.cpp

g++ -o Projection.os -c -O2 -pthread -w -O2 -pthread -fPIC -I. -I/usr/lib/jvm/java-6-sun/include -I/usr/lib/jvm/java-6-sun/include/linux  -I/home/juanjo/xmipp -I/home/juanjo/xmipp/libraries xmipp_Projection.cpp

g++ -o XmippJavaInterface.so -s -pthread -shared  -L/home/juanjo/xmipp/lib MetaData.os ImageDouble.os Projection.os -lXmippData

