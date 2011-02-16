

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XMIPP_EXTS {

    public final static String XMP = "xmp";
    public final static String DM3 = "dm3";
    public final static String SEL = "sel";
    public final static String PSD = "psd";
    public final static String MRCS = "mrcs";
    public final static String MRC = "mrc";
    public final static String IMG = "img";
    public final static String HED = "hed";
    public final static String INF = "inf";
    public final static String RAW = "raw";
    public final static String TIF = "tif";
    public final static String SPI = "spi";
    public final static String STK = "stk";
    public final static String VOL = "vol";
    public final static String SER = "ser";
    public final static String SPE = "spe";
    public final static String EXTS[] = {XMP, DM3, SEL, PSD, MRCS,
        MRC, IMG, HED, INF, RAW, TIF, SPI, STK, VOL, SER, SPE};

    public static boolean isSelFile(String filename) {
        return isType(filename, SEL);
    }

    public static boolean isSupportedType(String filename) {
        for (int i = 0; i < EXTS.length; i++) {
            if (isType(filename, EXTS[i])) {
                return true;
            }
        }

        return false;
    }

    private static boolean isType(String filename, String ext) {
        if (filename.toLowerCase().endsWith(ext.toLowerCase())) {
            return true;
        }

        return false;
    }
}
