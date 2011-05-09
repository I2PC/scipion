

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XMIPP_EXTS {

    public final static String DM3 = "dm3";
    public final static String HED = "hed";
    public final static String INF = "inf";
    public final static String IMG = "img";
    public final static String MRC = "mrc";
    public final static String MRCS = "mrcs";
    public final static String PSD = "psd";
    public final static String RAW = "raw";
    public final static String SEL = "sel";
    public final static String SER = "ser";
    public final static String SPE = "spe";
    public final static String SPI = "spi";
    public final static String STK = "stk";
    public final static String TIF = "tif";
    public final static String VOL = "vol";
    public final static String XMD = "xmd";
    public final static String XMP = "xmp";
    public final static String EXTS[] = {
        DM3, HED, INF, IMG, MRC, MRCS, PSD, RAW,
        SEL, SER, SPE, SPI, STK, TIF, VOL, XMD, XMP};
    public final static String SELEXTS[] = {SEL, XMD};

    public static boolean isSelFile(String filename) {
        return isType(filename, SELEXTS);
    }

    public static boolean isSupportedType(String filename) {
        return isType(filename, EXTS);
    }

    public static boolean isType(String filename, String exts[]) {
        for (int i = 0; i < exts.length; i++) {
            if (isType(filename, exts[i])) {
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
