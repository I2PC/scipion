echo "Compiling java sources..."
./extract_labels.py
javac HelloWorld.java xmipp/*.java
jar -Mcf XmippJavaInterface.jar xmipp
echo "Generating header files"
javah xmipp.ImageDouble xmipp.MetaData
echo "Compiling library"
g++ -o ImageDouble.os -c -O2 -pthread -w -O2 -pthread -fPIC -I. -I/usr/lib64/jvm/java-1.6.0/include -I/usr/lib64/jvm/java/include/linux  -I/home/josem/xmipp_current -I/home/josem/xmipp_current/libraries xmipp_ImageDouble.cpp
g++ -o MetaData.os -c -O2 -pthread -w -O2 -pthread -fPIC -I. -I/usr/lib64/jvm/java-1.6.0/include -I/usr/lib64/jvm/java/include/linux  -I/home/josem/xmipp_current -I/home/josem/xmipp_current/libraries xmipp_MetaData.cpp    
g++ -o XmippJavaInterface.so -s -pthread -shared  -L/home/josem/xmipp_current/lib MetaData.os ImageDouble.os -lXmippData

