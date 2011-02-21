package xmipp;

public class ImageDouble {

    private static final String PSD_EXTENSION = ".psd";
    //hold pointer to Image class in C++ space
    private long peer;
    private String filename;
    private int width, height, depth;
    private long nimages;

    // Initialize library.
    static {
        System.loadLibrary("XmippJavaInterface");
        storeIds();
    }

    //caching some ids
    private static native void storeIds();

    //functions to create images
    private native void create();

    // Destructor.
    private synchronized native void destroy();

    // Reading.
    private native void read(String filename, boolean readData, long nimage) throws Exception;

    private native void readPreview(String filename,
            int w, int h, int slice, long nimage) throws Exception;

    // Set data.
    public native void setData(int w, int h, int d, double data[]) throws Exception;

    // Writting.
    public native void write(String filename) throws Exception;

    // Data.
    public native double[] getData();

    public native double getPixel(int x, int y);

    public native void setPixel(int x, int y, double value);

    public native double getVoxel(int x, int y, int z);

    public native void setVoxel(int x, int y, int z, double value);

    public native void convertPSD();

    public native int getXsize();

    public native int getYsize();

    public native int getZsize();

    public native long getNsize();

    public native void printShape();

    //non-native functions
    //constructor
    public ImageDouble() {
        create();
    }

    public ImageDouble(String filename) throws Exception {
        this();
        this.filename = filename;
        read(filename);
    }

    //should be called by GarbageCollector before destroying
    @Override
    protected void finalize() {
        //System.out.println("Destroying image: " + this.filename);
        destroy();
    }

    public void readHeader(String filename) throws Exception {
        read(filename, false);
    }

    public void read(String filename) throws Exception {
        read(filename, true);
    }

    public void read(String filename, long nimage) throws Exception {
        read(filename, true, nimage);
    }

    private void read(String filename, boolean readData) throws Exception {
        read(filename, readData, 0);
    }

//    private void read(String filename, boolean readData, long nimage) throws Exception {
//        read_(filename, readData, nimage);
///*
//        storeData(filename);*/
//    }
    public void readPreview(String filename, int w, int h) throws Exception {
        readPreview(filename, w, h, 0, 0);
    }

    public void readPreview(String filename, int w, int h, int slice) throws Exception {
        readPreview(filename, w, h, slice, 0);
    }
//
//    public void readPreview(String filename, int w, int h, int slice, long nimage) throws Exception {
//        readPreview_(filename, w, h, slice, nimage);
///*
//        storeData(filename);*/
//    }

//    public native int[] getSize();
/*    private void storeData(String filename) {
    this.filename = filename;

    int dimensions[] = getDimensions_();

    width = dimensions[INDEX_WIDTH];
    height = dimensions[INDEX_HEIGHT];
    depth = dimensions[INDEX_DEPTH];
    nimages = getNImages_();
    }*/
    /*
    public int getWidth() {
    return width;
    }

    public int getHeight() {
    return height;
    }

    public int getDepth() {
    return depth;
    }

    public long getNimages() {
    return nimages;
    }*/
    public void setFilename(String filename) {
        this.filename = filename;
    }

    public boolean isPSD() {
        return filename.endsWith(PSD_EXTENSION);
    }
}
