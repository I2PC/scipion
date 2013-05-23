package xmipp.jni;

public class PickingClassifier
{
	public native void autopick(String micrograph);//generates micrograph pos file
	
	public native void correct(String micrograph);
	
	
	
}
