package xmipp;

public class ImageDouble {
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
	public native void read(String imgName);
	public native void write(String imgName);
	public native double[] getData();
	public native int[] getSize();
	public native void printShape();
	
	//non-native functions
	//constructor
	public ImageDouble()
	{
		create();
	}
	
	public ImageDouble(String imgName)
	{
		create();
		filename = imgName;
		read(imgName);
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