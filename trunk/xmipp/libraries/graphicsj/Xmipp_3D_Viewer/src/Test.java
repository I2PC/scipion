public class Test{
	static{
		System.loadLibrary("XmippDataJava");
	}

	public static void main(String args[]){
		Projection p=new Projection();
		p.reset(10,10);
	}
}
