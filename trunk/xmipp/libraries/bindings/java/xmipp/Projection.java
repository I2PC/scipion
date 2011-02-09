package xmipp;

public class Projection extends ImageDouble {

    private native void create();

    public native void reset(int h, int w);

    public native void setXmippOrigin();

    public static native void projectVolume(ImageDouble image, Projection projection,
            int height, int width,
            double rot, double tilt, double pshi);
}
