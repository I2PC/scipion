package xmipp;

public class MetaData {
	//hold pointer to Image class in C++ space
	private long peer;	
	private String filename;
	
	//caching some ids
	private static native void storeIds();
	
	//functions to create images
	private native void create();
	//destructor
	private synchronized native void destroy();
	//reading
	public native void read(String filename);
	public native void write(String filename);
	public native void print();
	public native boolean containsLabel(int label);
	//get values from metadata
	public native int getValueInt(int label, long objId);
	public native double getValueDouble(int label, long objId);
	public native String getValueString(int label, long objId);
	public native boolean getValueBoolean(int label, long objId);
	//set values 
	public native boolean setValueInt(int label, int value, long objId);
	public native boolean setValueDouble(int label, double value, long objId);
	public native boolean setValueString(int label, String value, long objId);
	public native boolean setValueBoolean(int label, boolean value, long objId);
	
	public native long[] findObjects();
	
	//non-native functions
	//constructor
	public MetaData()
	{
		create();
	}
	
	public MetaData(String filename)
	{
		create();
		this.filename = filename;
		read(filename);
	}
	

	//should be called by GarbageCollector before destroying
	protected void finalize() 
	{
		destroy();
	}
	
	static
	{
		System.loadLibrary("XmippJavaInterface");
		storeIds();
	}
}