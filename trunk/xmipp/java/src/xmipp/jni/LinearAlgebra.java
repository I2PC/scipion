package xmipp.jni;

public class LinearAlgebra {
	
		//Load the native library
		static
		{
			System.loadLibrary("XmippJNI");
			storeIds();
		}
		// hold pointer to Image class in C++ space
		private long peer;

		public native void solveLinearSystem();

		// caching some ids
		private static native void storeIds();
		
		public native void clear();

		// functions to create images
		private native void create();

		// destructor
		private synchronized native void destroy();

		// non-native functions
		// constructor
		public LinearAlgebra()
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
