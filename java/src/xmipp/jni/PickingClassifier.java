package xmipp.jni;

public class PickingClassifier
{

	// Initialize library.
	static
	{
		System.loadLibrary("XmippJNI");
	}

	public native void autopick(String micrograph);

	public native void correct(String micrograph);

	public native void train();

}
