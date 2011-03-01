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
    // @TODO Set nicer icons :)
    public static final String TABLE[][] = {
        {"dm3", "/resources/type_dm3.png"},
        {"hed", "/resources/type_image.png"},
        {"inf", "/resources/type_image.png"},
        {"img", "/resources/type_image.png"},
        {"mrc", "/resources/type_image.png"},
        {"mrcs", "/resources/type_image.png"},
        {"psd", "/resources/type_image.png"},
        {"raw", "/resources/type_image.png"},
        {"sel", "/resources/type_selfile.png"},
        {"ser", "/resources/type_tia.png"},
        {"spe", "/resources/type_spe.png"},
        {"spi", "/resources/type_image.png"},
        {"stk", "/resources/type_image.png"},
        {"tif", "/resources/type_image.png"},
        {"vol", "/resources/type_image.png"},
        {"xmp", "/resources/type_image.png"},};

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
