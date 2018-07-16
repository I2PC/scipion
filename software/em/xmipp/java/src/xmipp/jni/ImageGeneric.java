package xmipp.jni;

import java.io.File;

public class ImageGeneric {

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
    public final static int ULong = 7;               // Unsigned integer (4 or 8 byte; depending on system)
    public final static int Long = 8;               // Signed integer (4 or 8 byte; depending on system)
    public final static int Float = 9;              // Floating point (4-byte)
    public final static int Double = 10;             // Double precision floating point (8-byte)
    public final static int ComplexShort = 11;      // Complex two-byte integer (4-byte)
    public final static int ComplexInt = 12;        // Complex integer (8-byte)
    public final static int ComplexFloat = 13;      // Complex floating point (8-byte)
    public final static int ComplexDouble = 14;     // Complex floating point (16-byte)
    public final static int Bool = 15;              // Boolean (1-byte?)
    public final static int UHalfByte = 16 ;         // 4 bits per pixel
    public final static int LastEntry = 17;         // This must be the last entry
    
    //AxisView constants for volume reslice
    public final static int Z_NEG = 0;
    public final static int Z_POS = 1; 
    public final static int Y_NEG = 2;    // Align -Y axis to Z axis, rotating 90 degrees around X axis");
    public final static int Y_POS = 3; // Align Y axis to Z axis, rotating -90 degrees around X axis");
    public final static int X_NEG = 4;   // Align -X axis to Z axis, rotating -90 degrees around Y axis");
    public final static int X_POS = 5;  // Align X axis to Z axis, rotating 90 degrees around Y axis");

    public final static int VIEWS[] = {Z_NEG, Y_NEG, X_NEG, Y_POS, X_POS };

    
    
    // Associated filename.
    private String filename;
    private boolean useLogarithm = true;   // To convert PSD images.
    // pointer to Image class in C++ space. Needed by native library.
    long peer;

    // Initialize library.
    static {
        System.loadLibrary("XmippJNI");
        //storeIds();
    }

    // Catching some ids
    //private static native void storeIds();
    // Functions to create objects.
    private native void create();

    // Destructor.
    public synchronized native void destroy();

    //non-native functions
    //constructor
    public ImageGeneric() throws Exception {
        create();
    }

    public ImageGeneric(int datatype) throws Exception {
        this();

        setDataType(datatype);
    }

    public ImageGeneric(String filename) throws Exception {
        this();
        setFilename(filename);
        
    }
    
    /** Set a new filename and reload the header */
    public void setFilename(String filename) throws Exception {
    	this.filename = filename;
    	readHeader(filename);//Filename.getFilename(filename));
    }

    public void resize(int w, int h) throws Exception {
        resize(w, h, 1);
    }

    public void resize(int w, int h, int d) throws Exception {
        resize(w, h, d, 1);
    }

    public native void resize(int h, int w, int d, long n) throws Exception;

    public String getFilename() {
        return filename;
    }

    public native int getXDim() throws Exception;

    public native int getYDim() throws Exception;

    public native int getZDim() throws Exception;

    public native long getNDim() throws Exception;

    public native int getDataType() throws Exception;

    public boolean getUseLogarithm() {
        return useLogarithm;
    }

    public void setUseLogarithm(boolean useLogarithm) {
        this.useLogarithm = useLogarithm;
    }

    // Header reader.
    private native void readHeader(String filename) throws Exception;

    // Data readers.
    public void read(int slice) throws Exception {
        read(slice, FIRST_IMAGE);
    }

    public void read(int width, int height, int slice) throws Exception {
        read(width, height, slice, FIRST_IMAGE);
    }

    public void read(int slice, long image) throws Exception {
        read(getXDim(), getYDim(), slice, image);
    }

    public void read(long image) throws Exception {
        read(ALL_SLICES, image);
    }

    public void read(int width, int height, long image) throws Exception {
        read(width, height, ALL_SLICES, image);
    }

    public void read(int width, int height, int slice, long image) throws Exception {
        read(filename, width, height, slice, image, true);
        // At this moment we don't know why
//        read(Filename.getFilename(filename), width, height, slice, image, true);
    }
    
    public void read(String filename, boolean map) throws Exception {
    	setFilename(filename);
    	read(filename, getXDim(), getYDim(), MID_SLICE, ALL_IMAGES, map);//image in filename will be used    	
    }
    
    public void read(String filename, int width, int height) throws Exception {
        this.filename = filename;
    	read(filename, width, height, MID_SLICE, ALL_IMAGES, true);//image in filename will be used
    }
    
    public void read(int width, int height) throws Exception {
        read(filename, width, height, MID_SLICE, ALL_IMAGES, true);//image in filename will be used
    }
    
    private native void read(String filename, int width, int height, int slice, long nimage, boolean map) throws Exception;

    public void readApplyGeo(String filename, MetaData metadata, long id) throws Exception {
        readApplyGeo(filename, metadata, id, true);
    }
    
    public void readApplyGeo(String filename, MetaData metadata, long id, boolean wrap) throws Exception {
        //ImageGeneric image = new ImageGeneric(filename);
        readHeader(filename);
        readApplyGeo(filename, metadata, id, getXDim(), getYDim(), wrap);
    }
    
    public void readApplyGeo(String filename, MetaData metadata, long id, int w, int h) throws Exception {
        this.filename = filename;
        readApplyGeo_(filename, metadata, id, w, h, true);
    }
    
    public void readApplyGeo(String filename, MetaData metadata, long id, int w, int h, boolean wrap) throws Exception {
        this.filename = filename;
        readApplyGeo_(filename, metadata, id, w, h, wrap);
    }
    
    public void readApplyGeo(MetaData metadata, long id, int w, int h) throws Exception {
        readApplyGeo_(filename, metadata, id, w, h, true);
    }

    public void readApplyGeo(MetaData metadata, long id, int w, int h, boolean wrap) throws Exception {
        readApplyGeo_(filename, metadata, id, w, h, wrap);
    }
    
    private native void readApplyGeo_(String filename, MetaData metadata, long id, int w, int h, boolean wrap) throws Exception;

    // Getters for data arrays.
    public native byte[] getArrayByte(long select_image, int slice) throws Exception;

    public native short[] getArrayShort(long select_image, int slice) throws Exception;

    public native float[] getArrayFloat(long select_image, int slice) throws Exception;    

    // Setters for data arrays.
    public native void setArrayByte(byte data[], long select_image, int slice) throws Exception;

    public native void setArrayShort(short data[], long select_image, int slice) throws Exception;

    public native void setArrayFloat(float data[], long select_image, int slice) throws Exception;

    //Some others image generic utilities
    public final native void setDataType(int dataType) throws Exception;

    public final native void convert2Datatype(int dataType) throws Exception;

    public native void mapFile2Write(int xDim, int yDim, int zDim, String filename, long nimage) throws Exception;

    public native void write(String filename) throws Exception;

    public native void printShape() throws Exception;

    public native double[] getStatistics() throws Exception;

    public native void setXmippOrigin() throws Exception;

    public native void convertPSD(boolean useLogarithm) throws Exception;
   
    public native void generatePSDCTF(MetaData md) throws Exception;
    
    public native void generateImageWithTwoCTFs(MetaData md1, MetaData md2, int xdim) throws Exception;
   
    public native void getReslice(ImageGeneric imgOut, int view) throws Exception;
    
    public native void reslice(int view) throws Exception;
    
    
    public native void getPreview(ImageGeneric imgOut, int xdim, int ydim, 
    		int select_slice, long select_image) throws Exception;
    

    public synchronized native double[] alignImage(ImageGeneric img) throws Exception;
    
    public synchronized native void getRadialAvg(ImageGeneric radialimg) throws Exception;
    
    public native void applyAlignment(ImageGeneric img, int index, double rot, double tilt, double psi) throws Exception;

    
    public native void removeAlignment(ImageGeneric img, int index, double rot, double tilt, double psi) throws Exception;
    
    public native Particle bestShift(ImageGeneric img) throws Exception;


    //Check if two images have same values to some accuracy
    public native boolean equal(ImageGeneric img, double accuracy) throws Exception;

    //subtract two imageGeneric. calling image minuend
    public native void subtract(ImageGeneric imgSubtrahend, ImageGeneric imgResult) throws Exception;

    //Smooth image generic by color dithering
    public native void smooth(ImageGeneric imgResult) throws Exception;
    
    public boolean isPSD() {
        return Filename.isPSD(filename);
    }

    public boolean isStack() throws Exception {
        return getNDim() > 1;
    }

    public boolean isVolume() throws Exception {
        return getZDim() > 1;
    }

    public boolean isStackOrVolume() throws Exception {
        return isStack() || isVolume();
    }

    public boolean isSingleImage() throws Exception {
        return !isStackOrVolume();
    }

    // Should be called by GarbageCollector before destroying
    @Override
    @SuppressWarnings("FinalizeDeclaration")
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }

	public static boolean exists(String imagepath)
	{
		if(imagepath == null || imagepath.isEmpty())
			return false;
		String prefix = Filename.getPrefix(imagepath);
		if(prefix == null)
			return new File(imagepath).exists();
		try
		{
		
			int index = Integer.parseInt(prefix);
			String file = Filename.removePrefix(imagepath);
			ImageGeneric ig = new ImageGeneric(file);
			ig.readHeader(file);
			if(index < 0 || index > ig.getNDim())
				return false;
			return true;
						
		}
		catch (Exception e)
		{
			return false;
		}
	}
	
	public native void applyGeo(double shiftX, double shiftY, double anglepsi, boolean flip, boolean wrap, double scaleFactor);
        
    public void applyGeo(double shiftX, double shiftY, double anglepsi, boolean flip, boolean wrap)
    {
    	applyGeo(shiftX, shiftY, anglepsi, flip, wrap, 1.);
    }    
    /**
     * Same as previous function, but reading the geometry information
     * from the transformation matrix. The matrix will be passed to the 
     * binding in its string representation read from the sqlite file.
     */
    
    public void applyGeoMatrix(String matrixString, boolean wrap)
    {
    	applyGeoMatrix(matrixString, 1, wrap);
    }
    
    
    public native void applyGeoMatrix(String matrixString, double scaleFactor, boolean wrap);
    

    
}
