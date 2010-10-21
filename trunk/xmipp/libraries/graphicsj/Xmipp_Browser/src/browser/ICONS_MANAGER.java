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
    public static int CACHE_SIZE_BROWSER = 20;
    public static int CACHE_SIZE_TABLE = 50;
    public final static String PATH_ICONS = "/resources/";
    public final static String PATH_FILTERING_ALERT = PATH_ICONS + "filtering.png";
    public final static String PATH_PARENT_DIRECTORY = PATH_ICONS + "parent_folder.png";
    public final static String PATH_REFRESH_DIRECTORY = PATH_ICONS + "refresh.png";
    public final static String PATH_DIRECTORY_FILE_TYPE = PATH_ICONS + "folder.png";
    public final static String PATH_UNKNOWN_FILE_TYPE = PATH_ICONS + "ij-unknown.png";
    public final static String PATH_IMAGEJ_IMAGE_FILE_TYPE = PATH_ICONS + "ij-image.png";
    public final static String PATH_IMAGEJ_FILE_TYPE = PATH_ICONS + "ij-file.png";
    public final static String PATH_IMAGEJ_XMIPP_SELFILE_TYPE = PATH_ICONS + "selfile.png";
    public final static String PATH_MISSING_ITEM = PATH_ICONS + "missing.png";
    public final static String PATH_DISABLED_ITEM = PATH_ICONS + "disabled.png";
    public final static String PATH_CAPTURE_WINDOW = PATH_ICONS + "capture.png";
    public static ImageIcon FILTERING_ALERT;
    public static ImageIcon PARENT_DIRECTORY;
    public static ImageIcon REFRESH_DIRECTORY;
    public static ImageIcon DIRECTORY_FILE_TYPE;
    public static ImageIcon IMAGEJ_IMAGE_FILE_TYPE;
    public static ImageIcon IMAGEJ_FILE_TYPE;
    public static ImageIcon XMIPP_SELFILE_TYPE;
    public static ImageIcon UNKNOWN_FILE_TYPE;
    public static ImageIcon MISSING_ITEM;
    public static ImageIcon DISABLED_ITEM;
    public static ImageIcon CAPTURE_WINDOW;

    public ICONS_MANAGER() {
        FILTERING_ALERT = new ImageIcon(getClass().getResource(PATH_FILTERING_ALERT));
        PARENT_DIRECTORY = new ImageIcon(getClass().getResource(PATH_PARENT_DIRECTORY));
        REFRESH_DIRECTORY = new ImageIcon(getClass().getResource(PATH_REFRESH_DIRECTORY));
        DIRECTORY_FILE_TYPE = new ImageIcon(getClass().getResource(PATH_DIRECTORY_FILE_TYPE));
        IMAGEJ_IMAGE_FILE_TYPE = new ImageIcon(getClass().getResource(PATH_IMAGEJ_IMAGE_FILE_TYPE));
        IMAGEJ_FILE_TYPE = new ImageIcon(getClass().getResource(PATH_IMAGEJ_FILE_TYPE));
        XMIPP_SELFILE_TYPE = new ImageIcon(getClass().getResource(PATH_IMAGEJ_XMIPP_SELFILE_TYPE));
        UNKNOWN_FILE_TYPE = new ImageIcon(getClass().getResource(PATH_UNKNOWN_FILE_TYPE));
        MISSING_ITEM = new ImageIcon(getClass().getResource(PATH_MISSING_ITEM));
        DISABLED_ITEM = new ImageIcon(getClass().getResource(PATH_DISABLED_ITEM));
        CAPTURE_WINDOW = new ImageIcon(getClass().getResource(PATH_CAPTURE_WINDOW));
    }
}
