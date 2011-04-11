package xmipp;
/**
 * Protocol for integrating native C++ code - @see ImageDouble.java
 */
public class MetaData {

    //hold pointer to Image class in C++ space
    private long peer;
    private String filename;

    static {
        System.loadLibrary("XmippJavaInterface");
        storeIds();
    }

    //caching some ids
    private static native void storeIds();

    //functions to create images
    private native void create();

    //destructor
    private synchronized native void destroy();

    //reading
    public native void read(String filename);

    public native int size();

    public native void write(String filename);

    public native void print();

    public native boolean containsLabel(int label);

    //get values from metadata
    public native int getValueInt(int label, long objId);

    public native double getValueDouble(int label, long objId);

    public native String getValueString(int label, long objId);

    public native boolean getValueBoolean(int label, long objId);

    //set values
    public native boolean setValueInt(int label, int value, long objId);

    public native boolean setValueDouble(int label, double value, long objId);

    public native boolean setValueString(int label, String value, long objId);

    public native boolean setValueBoolean(int label, boolean value, long objId);

    public native long[] findObjects();

    public native long addObject();

    public native void addLabel(int label);

    public native String getFilename();

    //non-native functions
    //constructor
    public MetaData() {
        create();
    }

    public MetaData(String filename) {
        create();
        this.filename = filename;
        read(filename);
    }

    // Should be called by GarbageCollector before destroying
    @Override
    protected void finalize() throws Throwable {
	super.finalize();
        destroy();
    }
}
