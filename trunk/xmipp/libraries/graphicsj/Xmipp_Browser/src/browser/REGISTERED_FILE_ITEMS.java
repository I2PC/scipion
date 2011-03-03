/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser;

import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class REGISTERED_FILE_ITEMS {

    public static final int INDEX_EXT = 0;
    public static final int INDEX_ICON = 1;
    public static final String TABLE[][] = {
        {"jpeg", "/resources/type_image.png"},
        {"jpg", "/resources/type_image.png"},
        {"gif", "/resources/type_image.png"},
        {"png", "/resources/type_image.png"},
        {"dm3", "/resources/type_dm3.png"},
        {"hed", "/resources/type_hed.png"},
        {"inf", "/resources/type_inf.png"},
        {"img", "/resources/type_img.png"},
        {"mrc", "/resources/type_mrc.png"},
        {"mrcs", "/resources/type_mrcs.png"},
        {"psd", "/resources/type_psd.png"},
        {"raw", "/resources/type_raw.png"},
        {"sel", "/resources/type_sel.png"},
        {"ser", "/resources/type_ser.png"},
        {"spe", "/resources/type_spe.png"},
        {"spi", "/resources/type_spi.png"},
        {"stk", "/resources/type_stk.png"},
        {"tif", "/resources/type_tif.png"},
        {"vol", "/resources/type_vol.png"},
        {"xmd", "/resources/type_xmd.png"},
        {"xmp", "/resources/type_xmp.png"},};

    public static String getExtension(int index) {
        return TABLE[index][INDEX_EXT];
    }

    public static String getIcon(int index) {
        return TABLE[index][INDEX_ICON];
    }

    public static int length() {
        return TABLE.length;
    }

    public static String getFileExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (f.isDirectory()) {
            ext = null;
        } else if (i > 0 && i < s.length() - 1) {
            ext = s.substring(i + 1).toLowerCase();
        }

        return ext;
    }

    public static int getIndexForRegisteredItem(File file) {
        return getIndexForRegisteredItem(getFileExtension(file));
    }

    public static int getIndexForRegisteredItem(String ext) {
        if (ext != null) {
            for (int i = 0; i < TABLE.length; i++) {
                if (ext.compareTo(getExtension(i)) == 0) {
                    return i;
                }
            }
        }

        return -1;
    }
}
