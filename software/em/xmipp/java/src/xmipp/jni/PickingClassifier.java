package xmipp.jni;


public class PickingClassifier extends Classifier
{

    // pointer to AutoParticlePicking2 class in C++ space. Needed by native library.
    long peer;
    
	// Initialize library.
	static
	{
		System.loadLibrary("XmippJNI");
	}
	
	public PickingClassifier(int particlesize, String output, MDRow[] mics) throws Exception {
        create(particlesize, output, mics);
    }
	
	private native void create(int particle_size, String output, MDRow[] mics);
	
	public native void destroy();

	public synchronized native Particle[] autopick(String micrograph, int percent);
	
	public synchronized native void correct(MDRow[] manualRows, MDRow[] automaticRows);

	public synchronized native void train(MDRow[] micrographs, int x, int y, int width, int height);
	
	public synchronized native void setSize(int psize);
        
    public native int getParticlesThreshold();
	
	
	
    // Should be called by GarbageCollector before destroying
    @Override
    @SuppressWarnings("FinalizeDeclaration")
    protected void finalize() throws Throwable {
        super.finalize();
        destroy();
    }
    
    @Override
	public boolean needsTraining()
	{
		return true;
	}

	@Override
	public int getTrainingParticlesMinimum()
	{
		// TODO Auto-generated method stub
		return 15;
	}

}
