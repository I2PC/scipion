package xmipp.jni;

public class ProgTomographAlignment
{
	//Load the native library
	static
	{
		System.loadLibrary("XmippJNI");
		storeIds();
	}
	// hold pointer to Image class in C++ space
	private long peer;
    // Catching some ids
    private static native void storeIds();
    // Functions to create objects.
    protected native void create();
    
	public native void setInputFilename(String name);
	public native void setRoot(String name);
	public native void produceSideInfo();
	public native void run();
	public native void writeTransformations(String name);
	public native int getIteration();

	// destructor 
	private synchronized native void destroy();

	// non-native functions
	// constructor
	public ProgTomographAlignment()
	{
		create();
	}

	// Should be called by GarbageCollector before destroying
	@Override
	protected void finalize() throws Throwable
	{
		super.finalize();
		destroy();
	}
	
	
	
}
