package xmipp.jni;

import java.io.File;

/**
 * Protocol for integrating native C++ code - @see ImageDouble.java
 */
public class MetaData {
	/** Enum values with Labels possible types */
	public static final int LABEL_NOTYPE = -1;
	public static final int LABEL_INT = 0;
    public static final int LABEL_BOOL = 1;
    public static final int LABEL_DOUBLE = 2;
    public static final int LABEL_FLOAT = 3;
    public static final int LABEL_STRING = 4;
    public static final int LABEL_VECTOR = 5;
    public static final int LABEL_LONG = 6;
    public static final int LABEL_VECTOR_LONG = 7;
//
//    // Fields whose content is a path. They will be "fixed" conveniently.
//    private final static int PATHS_FIELDS[] = {
//        MDLabel.MDL_ASSOCIATED_IMAGE1,
//        MDLabel.MDL_ASSOCIATED_IMAGE2,
//        MDLabel.MDL_ASSOCIATED_IMAGE3,
//        MDLabel.MDL_IMAGE,
//        MDLabel.MDL_PSD,
//        MDLabel.MDL_CTFMODEL
//    };
//
//    static {
//        // Sorts it to use binary search later.
//        // (It's executed just for the first time)
//        Arrays.sort(PATHS_FIELDS);
//    }

    private final static int GEOMETRY_LABELS[] = {
        MDLabel.MDL_FLIP,
        MDLabel.MDL_ANGLEPSI,
        //        MDLabel.MDL_ANGLEROT,
        //        MDLabel.MDL_ANGLETILT,
        MDLabel.MDL_SHIFTX,
        MDLabel.MDL_SHIFTY, //        MDLabel.MDL_SHIFTZ
    };
    
    public static final int MD_OVERWRITE = 0;
    public static final int MD_APPEND = 1;
    private String filename;
    //hold pointer to Image class in C++ space
    private long peer;

    //keep labels for avoid read all the time
    int[] activeLabels;
    
    
    static {
        System.loadLibrary("XmippJNI");
        //storeIds();
    }

    //caching some ids
    //private static native void storeIds();
    //functions to create images
    private native void create();

    //destructor
    private synchronized native void destroy();
    

    //reading
    public native void read_(String filename) throws Exception;

    public final void read(String filename) throws Exception {
        this.filename = filename;
        read_(filename);
        activeLabels = getActiveLabels();
    }

    public native int size() throws Exception;

    public native void setColumnFormat(boolean format) throws Exception;

    public native void write(String filename) throws Exception;

    public native void writeBlock(String filename) throws Exception;

    public native void print() throws Exception;

    public native boolean containsLabel(int label) throws Exception;

    public boolean containsGeometryInfo() throws Exception {
        for (int i = 0; i < GEOMETRY_LABELS.length; i++) {
            if (containsLabel(GEOMETRY_LABELS[i])) {
                return true;
            }
        }
        return false;
    }

    public static native String label2Str(int label) throws Exception;
    
    /** Same of before but handling the exception */
    public static String getLabelName(int label){
    	try {
    		return label2Str(label);
    	}
    	catch (Exception e){
    		return null;
    	}
    }

    public static native String[] getBlocksInMetaDataFile(String filename) throws Exception;

    public native int[] getActiveLabels() throws Exception;

    public static native int getLabelType(int label) throws Exception;
    
    public static Class getLabelClass(int label) throws Exception {
    	int type = getLabelType(label);
    	switch (type) {
    	case LABEL_INT:
    		return Integer.class;
    	case LABEL_BOOL:
    		return Boolean.class;
    	case LABEL_FLOAT:
    		return Float.class;
    	case LABEL_DOUBLE:
    		return Double.class;
    	case LABEL_LONG:
    		return Long.class;
    	case LABEL_STRING:
    	case LABEL_VECTOR:
    	case LABEL_VECTOR_LONG:
    		return String.class;   	
    		
    	}
    	return null;
    }

    public static native boolean isTextFile(int label) throws Exception;

    public static native boolean isMetadata(int label) throws Exception;

    public static native boolean isCtfParam(int label) throws Exception;

    public static native boolean isImage(int label) throws Exception;

    public static native boolean isStack(int label) throws Exception;

    public static native boolean isMicrograph(int label) throws Exception;

    public static native boolean isPSD(int label) throws Exception;

    public native void makeAbsPath(int label) throws Exception;
    
    //get values from metadata
    public native int getValueInt(int label, long objId) throws Exception;

    public native long getValueLong(int label, long objId) throws Exception;

    public native double getValueDouble(int label, long objId) throws Exception;

    /** Return the value of some label as String */
    public native String getValueString(int label, long objId) throws Exception;

    public String getValueString(int label, long objId, boolean fixPaths) throws Exception {
        String value = getValueString(label, objId);

        // Try to fix paths.
        if (fixPaths && filename != null && isPathField(label)) {
            value = fixPath(value);
        }

        return value;
    }
    
    /** Return all values of the row as String[] 
     * @throws Exception */
    public String[] getRowValues(long objId) throws Exception{
    	String[] values = new String[activeLabels.length];
    	for (int i = 0; i < activeLabels.length; ++i)
    		values[i] = getValueString(activeLabels[i], objId);
    	return values;
    }

    public static boolean isPathField(int label) throws Exception {
        return isTextFile(label) || isMetadata(label)
                || isCtfParam(label) || isImage(label)
                || isStack(label) || isMicrograph(label) || isPSD(label);
    }

//    public static boolean isPathField(int label) {
//        return Arrays.binarySearch(PATHS_FIELDS, label) >= 0;
//    }
    public String fixPath(String value) {
        return Filename.fixPath(value, getBaseDir(), true);
    }

    public String fixPath(String value, String baseDir) {
        return Filename.fixPath(value, baseDir, false);
    }

    public native boolean getValueBoolean(int label, long objId) throws Exception;

    public String getFilename() {
        return filename != null ? Filename.getFilename(filename) : "";
    }
    
    public String getBlock() {
        return Filename.getPrefix(filename);
    }

    public String getPath() {
        return filename;
    }

    public String getBaseDir() {
        File f = new File(getFilename());
        f = new File(f.toURI().normalize().getPath());

        return f.getParent();
    }

    public native double[] getStatistics(boolean applyGeo) throws Exception;

    public native double[] getColumnValues(int label) throws Exception;

    //set functions conection with MetaData class in C++
    public native boolean setValueInt(int label, int value, long objId) throws Exception;

    public native boolean setValueDouble(int label, double value, long objId) throws Exception;

    public native boolean setValueString(int label, String value, long objId) throws Exception;

    public native boolean setValueBoolean(int label, boolean value, long objId) throws Exception;

    /** Obtain all the objects ids in the MetaData */
    public native long[] findObjects() throws Exception;
    
    /** Order the metadata by some label.
     * You can order ASCending or DESCending.
     */
    public native void sort(int sortLabel, boolean ascending) throws Exception;

    /** Import objects from other MetaData */
    public native void importObjects(MetaData from, long ids[]) throws Exception;

    /** Return the id of the first object */
    public native long firstObject() throws Exception;

    /** Add a new object entry and return new id */
    public native long addObject() throws Exception;

    /** Add new column to MetaData */
    public native void addLabel(int label) throws Exception;

    /** Get the average and std images, result is left on input image */
    public native void getStatsImages(ImageGeneric imageAvg, ImageGeneric imageStd
    		, boolean applyGeo) throws Exception;
    
    public native void getPCAbasis(ImageGeneric basis) throws Exception;

    public native void computeFourierStatistics(MetaData mdIn) throws Exception;

    public native void enableDebug() throws Exception;

    //non-native functions
    //constructor
    public MetaData() {
        create();
    }

    public MetaData(String filename) throws Exception {
        create();
        read(filename);
    }

    // Should be called by GarbageCollector before destroying
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }

    public native void readPlain(String file, String columns) throws Exception;
}
