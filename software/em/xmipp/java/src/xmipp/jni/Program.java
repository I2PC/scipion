package xmipp.jni;

/**
 * Protocol for integrating native C++ code - @see ImageDouble.java
 */
public class Program {
    public static native int runByName(String progName, String args) throws Exception;
}
