package xmipp.jni;

public class CTFDescription {

    public static final int BACKGROUND_NOISE = 0;
    public static final int ENVELOPE = 1;
    public static final int PSD = 2;
    public static final int CTF = 3;
    public double[][] profiles, avgprofiles;
    public double FMAX; // Calculated when a file is loaded.
    public double downsampling;
    // pointer to Image class in C++ space. Needed by native library.
    long peer;

    // Initialize library.
    static {
        System.loadLibrary("XmippJNI");
        //storeIds();
    }

    // Catching some ids
    //private static native void storeIds();

    // Functions to create objects: Internally calls to "Produce_side_info()"
    protected native void create();

    // Destructor.
    protected synchronized native void destroy();

    // Reading.
    public native void read_(String filename) throws Exception;

    public void read(String filename) throws Exception {
        read_(filename);
        MetaData md = new MetaData(filename);
        if (md.containsLabel(MDLabel.MDL_CTF_DOWNSAMPLE_PERFORMED))
        	downsampling = md.getValueDouble(MDLabel.MDL_CTF_DOWNSAMPLE_PERFORMED, md.firstObject());
        else
        	downsampling = 1.0;
        md.destroy();
        
        FMAX = getFMAX()/downsampling;
    }

    // Data.
    private native double getFMAX();

    private native double[][] CTFProfile(double angle, double FMAX, int samples);

    public void CTFProfile(double angle, int samples) {
        profiles = CTFProfile(angle, FMAX, samples);
    }

    private native double[][] CTFAverageProfile(double FMAX, int samples);

    public void CTFAverageProfiles(int samples) {
        avgprofiles = CTFAverageProfile(FMAX, samples);
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
