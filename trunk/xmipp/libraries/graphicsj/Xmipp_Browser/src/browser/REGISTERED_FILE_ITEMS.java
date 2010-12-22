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
    public static final int INDEX_CLASS = 1;
    public static final int INDEX_ICON = 2;
    public static final String TABLE[][] = {
        {"img", "browser.imageitems.listitems.SpiderItem", "/resources/type_image.png"},
        {"xmp", "browser.imageitems.listitems.SpiderItem", "/resources/type_image.png"},
        {"vol", "browser.imageitems.listitems.SpiderItem", "/resources/type_image.png"},
        {"spi", "browser.imageitems.listitems.SpiderItem", "/resources/type_image.png"},
        {"sel", "browser.imageitems.listitems.SelFileItem", "/resources/type_selfile.png"},
        {"dm3", "browser.imageitems.listitems.DM3FileItem", "/resources/type_dm3.png"},
        {"spe", "browser.imageitems.listitems.SPEFileItem", "/resources/type_spe.png"},
        {"ser", "browser.imageitems.listitems.SERFileItem", "/resources/type_tia.png"}};

    public static String getExtension(int index) {
        return TABLE[index][INDEX_EXT];
    }

    public static String getClass(int index) {
        return TABLE[index][INDEX_CLASS];
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
