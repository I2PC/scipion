package xmipp.jni;

public class Projection {

    /*public static native void projectVolume(ImageGeneric volume, ImageGeneric projection,
            double rot, double tilt, double psi) throws Exception;*/

    public static native void projectVolume(ImageGeneric volume, ImageGeneric projection,
            double matrix[]) throws Exception;

    public static native double entropyOtsuSegmentation(ImageGeneric volume, double percentile, boolean binarize) throws Exception;

//    public static native double[] eulerDirection2Angles(double vector[]) throws Exception;
    public static native double[] eulerMatrix2Angles(double matrix[]) throws Exception;
}
