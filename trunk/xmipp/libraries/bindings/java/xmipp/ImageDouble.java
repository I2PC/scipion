package xmipp;

public class ImageDouble {

    public final static int FIRST_IMAGE = 1;
    public final static int FIRST_SLICE = 1;
    public final static int ALL_IMAGES = 0;
    public final static int ALL_SLICES = 0;
    public final static int MID_SLICE = -1;
    private final static String PSD_EXTENSION = ".psd";
    // pointer to Image class in C++ space. Needed by native library.
    private long peer;
    private String filename;

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

    // Reading.
    private native void read_image(String filename, boolean readData, long nimage) throws Exception;

    private native void readApplyGeo(MetaData metadata, long id) throws Exception;

    private native void read_preview(String filename,
            int w, int h, int slice, long nimage) throws Exception;

    // Set data.
    public native void setData(int w, int h, int d, double data[]) throws Exception;

    // Writting.
    public native void write(String filename) throws Exception;

    // Data.
    public native double[] getData();

    public native void convertPSD();

    public native int getXsize();

    public native int getYsize();

    public native int getZsize();

    public native long getNsize();

    public native void setXmippOrigin() throws Exception;

    public native void printShape();

    //non-native functions
    //constructor
    public ImageDouble() {
        create();
    }

    public ImageDouble(String filename) throws Exception {
        this();
        read(filename);
    }

    public ImageDouble(MetaData metadata, long id) throws Exception {
        this();
        read(metadata, id);
    }

    // Should be called by GarbageCollector before destroying
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }

    public void readHeader(String filename) throws Exception {
        readHeader(filename, ALL_IMAGES);
    }

    public void readHeader(String filename, long nimage) throws Exception {
        read(filename, false, nimage);
    }

    public void read(String filename) throws Exception {
        read(filename, ALL_IMAGES);
    }

    public void read(String filename, long nimage) throws Exception {
        read(filename, true, nimage);
    }

    private void read(String filename, boolean readData, long nimage) throws Exception {
        read_image(filename, readData, nimage);
        setFilename(filename);
    }

    private void read(MetaData metadata, long id) throws Exception {
        readApplyGeo(metadata, id);

        if (metadata.containsLabel(MDLabel.MDL_IMAGE)) {
            setFilename(metadata.getValueString(MDLabel.MDL_IMAGE, id));
        }
    }

    public void readPreview(String filename, int w, int h) throws Exception {
        readPreview(filename, w, h, MID_SLICE);
    }

    public void readPreview(String filename, int w, int h, int slice) throws Exception {
        readPreview(filename, w, h, slice, FIRST_IMAGE);
    }

    public void readPreview(String filename, int w, int h, int slice, long nimage) throws Exception {
        read_preview(filename, w, h, slice, nimage);
        setFilename(filename);
    }

    public void setFilename(String filename) {
        this.filename = filename;
    }

    public String getFilename() {
        return filename;
    }

    public boolean isPSD() {
        return filename.endsWith(PSD_EXTENSION);
    }
}
