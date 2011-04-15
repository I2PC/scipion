/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser;

import javax.swing.ImageIcon;

/**
 *
 * @author Juanjo Vega
 */
public class ICONS_MANAGER {

    public static int PREVIEW_WIDTH = 80;
    public static int PREVIEW_HEIGHT = 80;
    //public static int CACHE_SIZE_BROWSER = 20;
    //public static int CACHE_SIZE_TABLE = 50;
    public final static String PATH_ICONS = "/resources/";
    public final static String PATH_FILTERING_ALERT = PATH_ICONS + "filtering.png";
    public final static String PATH_PARENT_DIRECTORY = PATH_ICONS + "parent_folder.png";
    public final static String PATH_REFRESH_DIRECTORY = PATH_ICONS + "refresh.png";
    public final static String PATH_MISSING_ITEM = PATH_ICONS + "missing.png";
    public final static String PATH_CAPTURE_WINDOW = PATH_ICONS + "capture.png";
    public final static String PATH_FILE_TYPE_DIRECTORY = PATH_ICONS + "type_folder.png";
    public final static String PATH_FILE_TYPE_UNKNOWN = PATH_ICONS + "type_unknown.png";
    public static ImageIcon FILTERING_ALERT;
    public static ImageIcon PARENT_DIRECTORY;
    public static ImageIcon REFRESH_DIRECTORY;
    public static ImageIcon MISSING_ITEM;
    public static ImageIcon CAPTURE_WINDOW;
    public static ImageIcon FILE_TYPE_DIRECTORY;
    public static ImageIcon FILE_TYPE_UNKNOWN;
    public static ImageIcon FILE_TYPES_ICONS[] = new ImageIcon[REGISTERED_FILE_ITEMS.length()];

    static {
        new ICONS_MANAGER();    // Auto-load.
    }

    public ICONS_MANAGER() {
        FILTERING_ALERT = new ImageIcon(getClass().getResource(PATH_FILTERING_ALERT));
        PARENT_DIRECTORY = new ImageIcon(getClass().getResource(PATH_PARENT_DIRECTORY));
        REFRESH_DIRECTORY = new ImageIcon(getClass().getResource(PATH_REFRESH_DIRECTORY));
        MISSING_ITEM = new ImageIcon(getClass().getResource(PATH_MISSING_ITEM));
        CAPTURE_WINDOW = new ImageIcon(getClass().getResource(PATH_CAPTURE_WINDOW));

        FILE_TYPE_DIRECTORY = new ImageIcon(getClass().getResource(PATH_FILE_TYPE_DIRECTORY));
        FILE_TYPE_UNKNOWN = new ImageIcon(getClass().getResource(PATH_FILE_TYPE_UNKNOWN));

        // Loads icons for file types.
        for (int i = 0; i < FILE_TYPES_ICONS.length; i++) {
            FILE_TYPES_ICONS[i] = new ImageIcon(getClass().getResource(REGISTERED_FILE_ITEMS.getIcon(i)));
        }
    }
}
