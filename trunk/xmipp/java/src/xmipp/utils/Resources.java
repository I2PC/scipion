/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.utils;

import ij.ImagePlus;
import javax.swing.ImageIcon;

/**
 *
 * @author Juanjo Vega
 */
public class Resources {

    public static int DEFAULT_PREVIEW_WIDTH = 80;
    public static int DEFAULT_PREVIEW_HEIGHT = 80;
    
    //Path to resource folder    
    //FIXME: Now hard coded
    public final static String PATH_ICONS = "/home/josem/xmipp_current/resources/";
    //Icon names    
    public final static String MISSING = "missing.png";
    public final static String WAIT = "wait.gif";
    public final static String WAIT_MENU = "wait_menu.gif";
    public final static String NORMALIZE = "histogram.png";
    public final static String APPLY_GEO = "rotate.png";
    public final static String ADJUST_COLS = "resize.png";

    public static ImageIcon WAIT_ICON;
    public static ImageIcon WAIT_MENU_ICON;
    public static ImagePlus MISSING_ITEM;

    static {
        new Resources();    // Auto-load.
    }
    
    // Create an icon using the xmipp resource path
    public static ImageIcon getIcon(String name){
    	return new ImageIcon(PATH_ICONS + name);
    }

    public Resources() {
        WAIT_ICON = getIcon(WAIT);
        WAIT_MENU_ICON = getIcon(WAIT_MENU);

        // For missing items.
        ImageIcon missingIcon = getIcon(MISSING);
        MISSING_ITEM = new ImagePlus("X", missingIcon.getImage());
    }
}
