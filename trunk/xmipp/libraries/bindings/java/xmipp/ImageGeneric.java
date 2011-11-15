package xmipp;

public class ImageGeneric {
	
	//Constants to image selection indexes
    public final static int FIRST_IMAGE = 1;
    public final static int FIRST_SLICE = 1;
    public final static int ALL_IMAGES = 0;
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

    public String filename;
    //Dimensions and datatype
    public int xSize;
    public int ySize;
    public int zSize;
    public long nSize;
    public int dataType;
    
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
    	create();       
        this.filename = filename;
        readHeader(filename);
    }

    
    private native void read(String filename, boolean onlyHeader);
    
    public native void write(String filename);
    
    public void readHeader(String filename) throws Exception {
    	read(filename, true);
    }
    
    public void readData(String filename) throws Exception {
    	read(filename, false);
    }
    
    public byte[] getArrayByte(){
    	return getArrayByte(filename, xSize, ySize, zSize, nSize, dataType);
    }
    
    private static native byte[] getArrayByte(String filename, int x, int y, int z, long N, int datatype);
    
    public short[] getArrayShort(){
    	return getArrayShort(filename, xSize, ySize, zSize, nSize, dataType);
    }
    
    private static native short[] getArrayShort(String filename, int x, int y, int z, long N, int datatype);
    
    public float[] getArrayFloat(){
    	return getArrayFloat(filename, xSize, ySize, zSize, nSize, dataType);
    }
    
    private static native float[] getArrayFloat(String filename, int x, int y, int z, long N, int datatype);
    
    public void setArrayFloat(float[] data) {
    	System.out.println("image info: x: " + xSize + " y:" + ySize);
    	System.out.println("data:");
		for (int j = 0; j < ySize; j++) {
			if (j < 3 || ySize-j <=3)
			{
				System.out.print("Line: " + j + " --> ");
				for (int i = 0; i < xSize; i++) {
					if (i < 3 || xSize-i<=3)
					System.out.print(data[j * ySize + i] + " ");
					else if (i==3)
						System.out.print("... ");
				}
				System.out.println("");
			}
			else if (j==3)
				System.out.println("...");
		}
		System.out.println("calling native setArrayFloat");
    	setArrayFloat(xSize, ySize, zSize, nSize, dataType, data);
    }
    
    private native void setArrayFloat(int x, int y, int z, long N, int datatype, float [] data);    

    public native double[] getStatistics() throws Exception;
    
    // Should be called by GarbageCollector before destroying
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }
    
   
}
