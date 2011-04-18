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
    public native void read(String filename);

    public native int size();

    public native void write(String filename);

    public native void print();

    public native boolean containsLabel(int label);

    //get values from metadata
    public native int getValueInt(int label, long objId);

    public native double getValueDouble(int label, long objId);

    private native String getValueString_(int label, long objId);

    public String getValueString(int label, long objId){
    	String value = getValueString_(label, objId);

	// Try to fix paths.
	if(Arrays.binarySearch(PATHS_FIELDS, label)>=0){
    		if(!value.startsWith(File.separator)){
    			value = fixPath(getBaseDir(), value); 
    		}
    	}

    	return value;
    }

    public native boolean getValueBoolean(int label, long objId);

    public native String getFilename();
    
    public String getBaseDir(){
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

    public MetaData(String filename) {
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
    private static String fixPath(String workdir, String filename) {
        int index;
        String aux = new String(filename.toCharArray());

        do {
            index = aux.lastIndexOf(File.separator);
            if (index >= 0) {
                aux = filename.substring(0, index);

                if (workdir.endsWith(aux)) {
                    filename = filename.substring(index);
                    break;
                }
            }
        } while (index >= 0);

        return workdir + filename;
    }
}
