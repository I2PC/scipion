package xmipp;

public class Projection {

    //hold pointer to Image class in C++ space
/*    private long peer;
    private String filename;*/

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

    public native void setXmippOrigin() throws Exception;

    public static native void projectVolume(ImageDouble image, Projection projection,
            double rot, double tilt, double pshi) throws Exception;

    // Writting.
    public native void write(String filename) throws Exception;

    // Data.
    public native double[] getData();

    public native int getXsize();

    public native int getYsize();

    public native int getZsize();

    public native void printShape();

    //non-native functions
    //constructor
    public Projection() {
        create();
    }
}
