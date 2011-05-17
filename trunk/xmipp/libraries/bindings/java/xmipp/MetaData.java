package xmipp;

import java.io.File;
import java.util.Arrays;

/**
 * Protocol for integrating native C++ code - @see ImageDouble.java
 */
public class MetaData {

    public final static String SEPARATOR = "@";
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
    public native void read(String filename) throws Exception;

    public native int size();

    public native void write(String filename) throws Exception;

    public native void print();

    public native boolean containsLabel(int label);

    //get values from metadata
    public native int getValueInt(int label, long objId);

    public native double getValueDouble(int label, long objId);

    private native String getValueString_(int label, long objId);

    public String getValueString(int label, long objId) {
        String value = getValueString_(label, objId);

        // Try to fix paths.
        if (Arrays.binarySearch(PATHS_FIELDS, label) >= 0) {
            value = fixPath(getBaseDir(), value);
        }

        return value;
    }

    public native boolean getValueBoolean(int label, long objId);

    public native String getFilename();

    public String getBaseDir() {
        File f = new File(getFilename());
        f = new File(f.getAbsolutePath());

        return f.getParent();
    }

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

    // Auxiliary methods.
    public static String fixPath(String workdir, String filename) {
        String fixed = filename;

        if (!filename.startsWith(File.separator)) { // Absolute path?
            String name = getFilename(filename);
            String strimage = "";

            if (filename.contains(SEPARATOR)) { // Has #image?
                int image = getNimage(filename);
                strimage = image + SEPARATOR;
            }

            if (!name.startsWith(File.separator)) { // In 'image@name', is name absolute?
                String aux = getAbsPath(workdir, name);
                fixed = strimage + aux;
            }
        }

        return fixed;
    }

    public static int getNimage(String filename) {
        int nimage = ImageDouble.ALL_IMAGES;

        if (filename.contains(SEPARATOR)) {
            String str = filename.split(SEPARATOR)[0];
            if (!str.isEmpty()) {
                nimage = Integer.valueOf(str);
            }
        }

        return nimage;
    }

    public static String getFilename(String filename) {
        if (filename.contains(SEPARATOR)) {
            return filename.split(SEPARATOR)[1];
        }

        return filename;
    }

    private static String getAbsPath(String baseDir, String filename) {
        String[] tokensDir = baseDir.split(File.separator);
        String[] tokensFile = filename.split(File.separator);

        int indexDir = tokensDir.length - 1;
        int indexFile = 0;

        while (indexFile < tokensFile.length && indexDir >= 0) {
            String dirToken = tokensDir[indexDir];
            String fileToken = tokensFile[indexFile];
            if (!dirToken.equals(fileToken)) {
                break;
            }
            indexDir--;
            indexFile++;
        }

        // Builds result path.
        String path = "";
        // Dir
        for (int i = 0; i < tokensDir.length; i++) {
            path += tokensDir[i] + File.separator;
        }

        // File
        for (int i = indexFile; i < tokensFile.length - 1; i++) {
            path += tokensFile[i] + File.separator;
        }
        path += tokensFile[tokensFile.length - 1];  // Last item (to avoid "/" at the end)

        return path;
    }
}
