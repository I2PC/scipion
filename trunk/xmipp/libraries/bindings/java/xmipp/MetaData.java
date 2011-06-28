package xmipp;

import java.io.File;
import java.util.Arrays;

/**
 * Protocol for integrating native C++ code - @see ImageDouble.java
 */
public class MetaData {

    // Fields whose content is a path. They will be "fixed" conveniently.
    private final static int PATHS_FIELDS[] = {
        MDLabel.MDL_ASSOCIATED_IMAGE1,
        MDLabel.MDL_ASSOCIATED_IMAGE2,
        MDLabel.MDL_ASSOCIATED_IMAGE3,
        MDLabel.MDL_IMAGE,
        MDLabel.MDL_PSD,
        MDLabel.MDL_CTFMODEL
    };

    static {
        // Sorts it to use binary search later.
        // (It's executed just for the first time)
        Arrays.sort(PATHS_FIELDS);
    }
    private String filename;
    //hold pointer to Image class in C++ space
    private long peer;

    static {
        System.loadLibrary("XmippJavaInterface");
        storeIds();
    }

    //caching some ids
    private static native void storeIds();

    //functions to create images
    private native void create();

    //destructor
    private synchronized native void destroy();

    //reading
    public native void read_(String filename) throws Exception;

    public void read(String filename) throws Exception {
        this.filename = filename;
        read_(filename);
    }

    public native int size();

    public native void write(String filename) throws Exception;

    public native void print();

    public native boolean containsLabel(int label);

    public static native String label2Str(int label);

    public native int[] getActiveLabels();

    public static native Class getLabelType(int label);

    //get values from metadata
    public native int getValueInt(int label, long objId);

    public native double getValueDouble(int label, long objId);

    public native String getValueString(int label, long objId);

    public String getValueString(int label, long objId, boolean fixPaths) {
        String value = getValueString(label, objId);

        // Try to fix paths.
        if (fixPaths && isPathField(label)) {
            value = fixPath(value);
        }

        return value;
    }

    public static boolean isPathField(int label) {
        return Arrays.binarySearch(PATHS_FIELDS, label) >= 0;
    }

    public String fixPath(String value) {
        return Filename.fixPath(getBaseDir(), value);
    }

    public native boolean getValueBoolean(int label, long objId);

    public String getFilename() {
        return filename;
    }

    public String getBaseDir() {
        File f = new File(getFilename());
        f = new File(f.getAbsolutePath());

        return f.getParent();
    }

    public native double[] getStatistics(boolean applyGeo);

    //set values
    public native boolean setValueInt(int label, int value, long objId);

    public native boolean setValueDouble(int label, double value, long objId);

    public native boolean setValueString(int label, String value, long objId);

    public native boolean setValueBoolean(int label, boolean value, long objId);

    public native long[] findObjects();

    public native long addObject();

    public native void addLabel(int label);

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
}
