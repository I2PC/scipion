package xmipp;

public class ImageGeneric_ {

    //Constants to image selection indexes
    public final static long FIRST_IMAGE = 1;
    public final static int FIRST_SLICE = 1;
    public final static long ALL_IMAGES = 0;
    public final static int ALL_SLICES = 0;
    public final static int MID_SLICE = -1;
    //Datatypes constants
    public final static int Default = -1;           // For writing purposes
    public final static int Unknown_Type = 0;       // Undefined data type
    public final static int UChar = 1;              // Unsigned character or byte type
    public final static int SChar = 2;              // Signed character (for CCP4)
    public final static int UShort = 3;             // Unsigned integer (2-byte)
    public final static int Short = 4;              // Signed integer (2-byte)
    public final static int UInt = 5;               // Unsigned integer (4-byte)
    public final static int Int = 6;                // Signed integer (4-byte)
    public final static int Long = 7;               // Signed integer (4 or 8 byte; depending on system)
    public final static int Float = 8;              // Floating point (4-byte)
    public final static int Double = 9;             // Double precision floating point (8-byte)
    public final static int ComplexShort = 10;      // Complex two-byte integer (4-byte)
    public final static int ComplexInt = 11;        // Complex integer (8-byte)
    public final static int ComplexFloat = 12;      // Complex floating point (8-byte)
    public final static int ComplexDouble = 13;     // Complex floating point (16-byte)
    public final static int Bool = 14;              // Boolean (1-byte?)
    public final static int LastEntry = 15;         // This must be the last entry
    // Associated filename.
    private String filename;
    // Dimensions and datatype
    private int xDim;
    private int yDim;
    private int zDim;
    private long nDim;
    private int dataType;
    private boolean useLogarithm = true;   // To convert PSD images.
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
    protected final native void create();

    // Destructor.
    protected synchronized native void destroy();

    //non-native functions
    //constructor
    public ImageGeneric_() {
        create();
    }

    public ImageGeneric_(String filename) throws Exception {
        create();

        this.filename = filename;
        readHeader(filename);
    }

    public String getFilename() {
        return filename;
    }

    public int getXDim() {
        return xDim;
    }

    public int getYDim() {
        return yDim;
    }

    public int getZDim() {
        return zDim;
    }

    public long getNDim() {
        return nDim;
    }

    public int getDataType() {
        return dataType;
    }

    public boolean getUseLogarithm() {
        return useLogarithm;
    }

    public void setUseLogarithm(boolean useLogarithm) {
        this.useLogarithm = useLogarithm;
    }

    // Header reader.
    private native void readHeader(String filename) throws Exception;

    // Data readers.
//    public void read() throws Exception {
//        read(ALL_SLICES, FIRST_IMAGE);
//    }
    public void read(int slice) throws Exception {
        read(slice, FIRST_IMAGE);
    }

    public void read(int width, int height, int slice) throws Exception {
        read(width, height, slice, FIRST_IMAGE);
    }

    public void read(int slice, long image) throws Exception {
        read(xDim, yDim, slice, image);
    }

    public void read(long image) throws Exception {
        read(ALL_SLICES, image);
    }

    public void read(int width, int height, long image) throws Exception {
        read(width, height, ALL_SLICES, image);
    }

    public void read(int width, int height, int slice, long image) throws Exception {
        read(filename, width, height, slice, image);
    }

    private native void read(String filename, int width, int height, int slice, long nimage) throws Exception;

    // Getters for data arrays.
    public byte[] getArrayByte(int slice) throws Exception {
        return getArrayByte(slice, dataType);
    }

    private native byte[] getArrayByte(int slice, int dataType) throws Exception;

    public short[] getArrayShort(int slice) throws Exception {
        return getArrayShort(slice, dataType);
    }

    private native short[] getArrayShort(int slice, int dataType) throws Exception;

    public float[] getArrayFloat(int slice) throws Exception {
        return getArrayFloat(slice, dataType);
    }

    private native float[] getArrayFloat(int slice, int dataType) throws Exception;
//
//    // TODO: Remove deprecated methods.
//    // @Deprecated
//    public byte[] getArrayByte() throws Exception {
//        return getArrayByte(filename, xDim, yDim, zDim, nDim, dataType);
//    }
//
//    // @Deprecated
//    private static native byte[] getArrayByte(String filename, int x, int y, int z, long N, int datatype) throws Exception;
//
//    // @Deprecated
//    public short[] getArrayShort() throws Exception {
//        return getArrayShort(filename, xDim, yDim, zDim, nDim, dataType);
//    }
//
//    // @Deprecated
//    private static native short[] getArrayShort(String filename, int x, int y, int z, long N, int datatype) throws Exception;

    // @Deprecated
//    public float[] getArrayFloat() throws Exception {
//        return getArrayFloat(filename, xSize, ySize, zSize, nSize, dataType);
//    }
    // @Deprecated
//    private static native float[] getArrayFloat(String filename, int x, int y, int z, long N, int datatype) throws Exception;
    // Writer.
    public native void write(String filename) throws Exception;

    // Setters for data arrays.
//    public void setArrayFloat(float[] data) throws Exception {
//        System.out.println("image info: x: " + xDim + " y:" + yDim);
//        System.out.println("data:");
//        for (int j = 0; j < yDim; j++) {
//            if (j < 3 || yDim - j <= 3) {
//                System.out.print("Line: " + j + " --> ");
//                for (int i = 0; i < xDim; i++) {
//                    if (i < 3 || xDim - i <= 3) {
//                        System.out.print(data[j * yDim + i] + " ");
//                    } else if (i == 3) {
//                        System.out.print("... ");
//                    }
//                }
//                System.out.println("");
//            } else if (j == 3) {
//                System.out.println("...");
//            }
//        }
//        System.out.println("calling native setArrayFloat");
//        setArrayFloat(xDim, yDim, zDim, nDim, dataType, data);
//    }
//    private native void setArrayFloat(int x, int y, int z, long N, int datatype, float data[]) throws Exception;
    public native void printShape() throws Exception;

    public native double[] getStatistics() throws Exception;

    public native void setXmippOrigin() throws Exception;

    public native void convertPSD(boolean useLogarithm);

    public boolean isPSD() {
        return Filename.isPSD(filename);
    }

    public boolean isStack() {
        return getNDim() > 1;
    }

    public boolean isVolume() {
        return getZDim() > 1;
    }

    public boolean isStackOrVolume() {
        return isStack() || isVolume();
    }

    public boolean isSingleImage() {
        return !isStackOrVolume();
    }

    // Should be called by GarbageCollector before destroying
    @Override
    @SuppressWarnings("FinalizeDeclaration")
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }
}
