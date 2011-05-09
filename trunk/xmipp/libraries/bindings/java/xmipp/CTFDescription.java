package xmipp;

public class CTFDescription {

    public static final int NOISE = 0;
    public static final int NOISE_ERROR2 = 1;
    public static final int NOISE_CTF2 = 2;
    public static final int NOISE_PURE = 3;
    public double[][] profiles;
    public double FMAX; // Calculated when a file is loaded.
    // pointer to Image class in C++ space. Needed by native library.
    private long peer;

    // Initialize library.
    static {
        System.loadLibrary("XmippJavaInterface");
        storeIds();
    }

    // Catching some ids
    private static native void storeIds();

    // Functions to create objects: Internally calls to "Produce_side_info()"
    protected native void create();

    // Destructor.
    protected synchronized native void destroy();

    // Reading.
    public native void read_(String filename) throws Exception;

    public void read(String filename) throws Exception {
        read_(filename);
        FMAX = getFMAX();
    }

    // Data.
    public native double getFMAX();

    private native double[][] CTFProfile(double angle, double FMAX, int samples);

    public void getCTFProfile(double angle, int samples) {
        profiles = CTFProfile(angle, FMAX, samples);
    }

    //non-native functions
    //constructor
    public CTFDescription() {
        create();
    }

    public CTFDescription(String filename) throws Exception {
        this();
        read(filename);
    }

    // Should be called by GarbageCollector before destroying
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }
}
