package xmipp.jni;


public class PickingClassifier
{

    // pointer to AutoParticlePicking2 class in C++ space. Needed by native library.
    long peer;
    
	// Initialize library.
	static
	{
		System.loadLibrary("XmippJNI");
	}
	
	public PickingClassifier(int particlesize, String output) throws Exception {
        create(particlesize, output);
    }
	
	private native void create(int particle_size, String output);
	
	public synchronized native void destroy();

	public native void autopick(String micrograph, MetaData outputmd, int percent);

	public native void correct(MetaData manualmd, MetaData automaticmd);

	public native void train(MetaData micrographs);
	
	public native void setSize(int psize);
	
    // Should be called by GarbageCollector before destroying
    @Override
    @SuppressWarnings("FinalizeDeclaration")
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }

}
