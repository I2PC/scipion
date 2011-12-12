package xmipp;

import java.io.File;

/**
 * Protocol for integrating native C++ code - @see ImageDouble.java
 */
public class MetaData {
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

    static {
        System.loadLibrary("XmippJavaInterface");
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
    }

    public native int size();

    public native void setColumnFormat(boolean format);

    public native void write(String filename) throws Exception;

    public native void writeBlock(String filename) throws Exception;

    public native void print();

    public native boolean containsLabel(int label);

    public boolean containsGeometryInfo() {
        for (int i = 0; i < GEOMETRY_LABELS.length; i++) {
            if (containsLabel(GEOMETRY_LABELS[i])) {
                return true;
            }
        }
        return false;
    }

    public static native String label2Str(int label);

    public static native String[] getBlocksInMetaDataFile(String filename) throws Exception;

    public native int[] getActiveLabels();

    public static native Class getLabelType(int label);

    public static native boolean isTextFile(int label);

    public static native boolean isMetadata(int label);

    public static native boolean isCtfParam(int label);

    public static native boolean isImage(int label);

    public static native boolean isStack(int label);

    public static native boolean isMicrograph(int label);

    public static native boolean isPSD(int label);

    //get values from metadata
    public native int getValueInt(int label, long objId);

    public native long getValueLong(int label, long objId);

    public native double getValueDouble(int label, long objId);

    public native String getValueString(int label, long objId);

    public String getValueString(int label, long objId, boolean fixPaths) {
        String value = getValueString(label, objId);

        // Try to fix paths.
        if (fixPaths && filename != null && isPathField(label)) {
            value = fixPath(value);
        }

        return value;
    }

    public static boolean isPathField(int label) {
        return isTextFile(label) || isMetadata(label)
                || isCtfParam(label) || isImage(label)
                || isStack(label) || isMicrograph(label) || isPSD(label);
    }

//    public static boolean isPathField(int label) {
//        return Arrays.binarySearch(PATHS_FIELDS, label) >= 0;
//    }
    public String fixPath(String value) {
        return Filename.fixPath(value, getBaseDir(),true);
    }
    
    public String fixPath(String value,String baseDir) {
        return Filename.fixPath(value, baseDir,false);
    }


    public native boolean getValueBoolean(int label, long objId);

    public String getFilename() {
        return Filename.getFilename(filename);
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

    public native double[] getStatistics(boolean applyGeo);

    public native double[] getColumnValues(int label) throws Exception;

    //set values
    public native boolean setValueInt(int label, int value, long objId);

    public native boolean setValueDouble(int label, double value, long objId);

    public native boolean setValueString(int label, String value, long objId);

    public native boolean setValueBoolean(int label, boolean value, long objId);

    public native long[] findObjects();

    public native void importObjects(MetaData from, long ids[]) throws Exception;

    public native long firstObject();

    public native long addObject();

    public native void addLabel(int label);

    public native void getPCAbasis(ImageGeneric basis) throws Exception;

    public native void computeFourierStatistics(String filename) throws Exception;

    public native void enableDebug();

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

    public native void readPlain(String file, String columns);
}
