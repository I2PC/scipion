package xmipp;

public class ImageGeneric {

    protected String filename;
    // pointer to Image class in C++ space. Needed by native library.
    private long peer;

    // Initialize library.
    static {
        System.loadLibrary("XmippJavaInterface");
        storeIds();
    }

    // Catching some ids
    private static native void storeIds();

    // Functions to create objects.
    protected native void create();

    // Destructor.
    protected synchronized native void destroy();

    //non-native functions
    //constructor
    public ImageGeneric() {
        create();
    }

    public ImageGeneric(String filename) throws Exception {
        this();

        this.filename = filename;
        read(filename);
    }

    public native void read(String filename) throws Exception;

    public native double[] getStatistics() throws Exception;

    // Should be called by GarbageCollector before destroying
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }
}
