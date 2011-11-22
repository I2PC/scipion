package xmipp;

/**
 * Protocol for integrating native C++ code
 * 
 * 1) define the Java method prototype here. For example, public native void write(String filename,int select_img, boolean isStack,
 *           int mode, boolean adjust) throws Exception;
 *          
 * 2) Call make to generate the C++ header (in this case, xmipp_ImageDouble.h)
 * 
 * 3) Search for the C++ method prototype in this C++ header (the signature probably will be rewritten to handle overloading and parameters) and
 *    copy the method prototype as is to xmipp_ImageDouble.cpp . This is important because if the method prototype
 *    differs, the JNI call will fail with a UnsatisfiedLinkError. Following the example, it would be 
 *    JNIEXPORT void JNICALL Java_xmipp_ImageDouble_write__Ljava_lang_String_2IZIZ(JNIEnv *, jobject, jstring, jint, jboolean, jint, jboolean);
 *    (quite a different name from the original "write")
 *    
 * 4) Add parameter names to this prototype in the cpp file (for clarity). Then code the method body
 * (you can reuse other methods as template, since most of the code is the same)
 * 
 * 5) Build the code. The result splits between the dynamic library (libXmippJavaInterface.so) and the Java classes
 * (XmippJavaInterface.jar) The library must be in your library path (LD_LIBRARY_PATH), The jar must be in
 * ImageJ plugins subdirectory
 *
 */
public class ImageDouble {

    public final static int FIRST_IMAGE = 1;
    public final static int FIRST_SLICE = 1;
    public final static int ALL_IMAGES = 0;
    public final static int ALL_SLICES = 0;
    public final static int MID_SLICE = -1;
    private String filename;
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

    // Reading.
    private native void read_image(String filename, boolean readData, long nimage) throws Exception;

    private native void read_preview(String filename,
            int w, int h, int slice, long nimage) throws Exception;

    private void readApplyGeo(String filename, MetaData metadata, long id) throws Exception {
        readApplyGeo(filename, metadata, id, getXsize(), getYsize());
    }

    private native void readApplyGeo(String filename, MetaData metadata, long id, int w, int h) throws Exception;

    // Set data.
    public void setData(int width, int height, int depth, double data[]) throws Exception {
        setData(width, height, depth, 1, data);
    }

    /**
     * Important: if you don't use some dimension (for example, depth) set it to 1
     */
    public native void setData(int width, int height, int depth, int numberOfSlices, double data[]) throws Exception;

    // Writing.
    public void write(String filename) throws Exception {
        write(filename, ALL_IMAGES, isStack(), ImageWriteMode.WRITE_OVERWRITE, CastWriteMode.CW_CAST);
    }

    public native void write(String filename, int select_img, boolean isStack, int mode, int castWriteMode) throws Exception;

    // Acess to data.
    public native double[] getData();

    public native double[] getData(long image, int slice);

    public native void convertPSD(boolean useLogarithm);

    public native int getXsize();

    public native int getYsize();

    public native int getZsize();

    public native long getNsize();

    public native void setXmippOrigin() throws Exception;

    public native void printShape();

    public native void fastEstimateEnhancedPSD(String filename, double downsampling, int w, int h) throws Exception;

    //non-native functions
    //constructor
    public ImageDouble() {
        create();
    }

    public ImageDouble(String filename) throws Exception {
        this();

        String name = Filename.getFilename(filename);
        long nimage = Filename.getNimage(filename);

        read(name, nimage);
    }

    public ImageDouble(String filename, long nimage) throws Exception {
        this();
        read(filename, nimage);
    }

    public ImageDouble(MetaData metadata, long id) throws Exception {
        this();

        String imgFilename = metadata.getValueString(MDLabel.MDL_IMAGE, id, true);
        read(imgFilename, metadata, id);
    }

    public ImageDouble(String filename, MetaData metadata, long id) throws Exception {
        this();
        read(filename, metadata, id);
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

    public void readStack(String filename) throws Exception {
        read(filename, ALL_IMAGES);
    }

    public void readSlice(String filename) throws Exception {
        // select_img/nimage is only used for reading all the images, otherwise
        // the slice is encoded in the filename
        read(filename, -10);
    }

    public void read(String filename, long nimage) throws Exception {
        read(filename, true, nimage);
    }

    private void read(String filename, boolean readData, long nimage) throws Exception {
        read_image(filename, readData, nimage);
        setFilename(filename);
    }

    private void read(String filename, MetaData metadata, long id) throws Exception {
        readApplyGeo(filename, metadata, id);

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

    public void readPreview(String filename, MetaData metadata, long id, int w, int h) throws Exception {
        readApplyGeo(filename, metadata, id, w, h);

        if (metadata.containsLabel(MDLabel.MDL_IMAGE)) {
            setFilename(metadata.getValueString(MDLabel.MDL_IMAGE, id));
        }
    }

    public void setFilename(String filename) {
        this.filename = filename;
    }

    public String getFilename() {
        return filename;
    }

    public boolean isPSD() {
        return Filename.isPSD(filename);
    }

    public boolean isStack() {
        return getNsize() > 1;
    }

    public boolean isVolume() {
        return getZsize() > 1;
    }

    public boolean isStackOrVolume() {
        return isStack() || isVolume();
    }

    public boolean isSingleImage() {
        return !isStackOrVolume();
    }
}
