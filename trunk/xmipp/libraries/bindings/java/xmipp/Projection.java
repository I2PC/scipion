package xmipp;

public class Projection {

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

    public native void reset(int h, int w) throws Exception;

    public static native void projectVolume(ImageGeneric_ image, Projection projection,
            double rot, double tilt, double pshi) throws Exception;

    public static native double entropyOtsuSegmentation(ImageGeneric_ volume, double percentile, boolean binarize) throws Exception;
    
    // Writting.
    public native void write(String filename) throws Exception;

    // Data.
    public native double[] getData();

    public native double[] getData(long nimage, int slice);

    public native int getXsize();

    public native int getYsize();

    public native int getZsize();

    public native void printShape();

    //non-native functions
    //constructor
    public Projection() {
        create();
    }

    // Should be called by GarbageCollector before destroying
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }
}
