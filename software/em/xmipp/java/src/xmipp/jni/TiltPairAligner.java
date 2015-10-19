package xmipp.jni;

public class TiltPairAligner
{
	//Load the native library
	static
	{
		System.loadLibrary("XmippJNI");
		storeIds();
	}
	// hold pointer to Image class in C++ space
	private long peer;

	public native void addParticleToAligner(int x1, int y1, int x2, int y2);

	public native Particle getTiltedParticle(int x1, int y1);
	
	public native Particle getUntiltedParticle(int x1, int y1);
	
	public native double[] computeAngles();

	// caching some ids
	private static native void storeIds();
	
	public native void clear();

	// functions to create images
	private native void create();

	// destructor
	private synchronized native void destroy();

	// non-native functions
	// constructor
	public TiltPairAligner()
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
