package xmipp.ij.plugins.maskstoolbar;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import javax.swing.ImageIcon;

import xmipp.jni.Filename;

/**
 *
 * @author Juanjo Vega
 */
@SuppressWarnings("ResultOfObjectAllocationIgnored")
public class ICONS {

    public final static String PATH_ICONS = Filename.getXmippPath("resources/");
    public final static String PATH_TOOL_RECTANGLE = PATH_ICONS + "rectangle.png";
    public final static String PATH_TOOL_ROUNDRECTANGLE = PATH_ICONS + "roundrectangle.png";
    public final static String PATH_TOOL_ELLIPSE = PATH_ICONS + "ellipse.png";
    public final static String PATH_TOOL_OVAL = PATH_ICONS + "oval.png";
    public final static String PATH_TOOL_POLYGON = PATH_ICONS + "polygon.png";
    public final static String PATH_TOOL_FREEHAND = PATH_ICONS + "freehand.png";
    public final static String PATH_TOOL_BRUSH = PATH_ICONS + "brush.png";
    public final static String PATH_ACTION_CREATEMASK = PATH_ICONS + "create_mask.png";
    public final static String PATH_ACTION_CREATESELECTION = PATH_ICONS + "create_selection.png";
    public final static String PATH_ACTION_INVERTSELECTION = PATH_ICONS + "invert_selection.png";
    public final static String PATH_ACTION_INVERTMASK = PATH_ICONS + "invert_mask.png";
    public final static String PATH_ACTION_SPECIFYSELECTION = PATH_ICONS + "specify_selection.png";
    public final static String PATH_ACTION_BRUSHSIZE = PATH_ICONS + "brush_size.png";
    public final static String PATH_MASK_LOCKED = PATH_ICONS + "locked.png";
    public final static String PATH_MASK_UNLOCKED = PATH_ICONS + "unlocked.png";
    public static ImageIcon TOOL_RECTANGLE;
    public static ImageIcon TOOL_ROUNDRECTANGLE;
    public static ImageIcon TOOL_ELLIPSE;
    public static ImageIcon TOOL_OVAL;
    public static ImageIcon TOOL_POLYGON;
    public static ImageIcon TOOL_FREEHAND;
    public static ImageIcon TOOL_BRUSH;
    public static ImageIcon ACTION_CREATEMASK;
    public static ImageIcon ACTION_CREATESELECTION;
    public static ImageIcon ACTION_INVERTSELECTION;
    public static ImageIcon ACTION_INVERTMASK;
    public static ImageIcon ACTION_SPECIFYSELECTION;
    public static ImageIcon ACTION_BRUSHSIZE;
    public static ImageIcon MASK_LOCKED;
    public static ImageIcon MASK_UNLOCKED;

    static {
        new ICONS();    // Auto-load.
    }

    public ICONS() {
//        TOOL_RECTANGLE = new ImageIcon(getClass().getResource(PATH_TOOL_RECTANGLE));
//        TOOL_ROUNDRECTANGLE = new ImageIcon(getClass().getResource(PATH_TOOL_ROUNDRECTANGLE));
//        TOOL_ELLIPSE = new ImageIcon(getClass().getResource(PATH_TOOL_ELLIPSE));
//        TOOL_OVAL = new ImageIcon(getClass().getResource(PATH_TOOL_OVAL));
//        TOOL_POLYGON = new ImageIcon(getClass().getResource(PATH_TOOL_POLYGON));
//        TOOL_FREEHAND = new ImageIcon(getClass().getResource(PATH_TOOL_FREEHAND));
//        TOOL_BRUSH = new ImageIcon(getClass().getResource(PATH_TOOL_BRUSH));
//        ACTION_CREATEMASK = new ImageIcon(getClass().getResource(PATH_ACTION_CREATEMASK));
//        ACTION_CREATESELECTION = new ImageIcon(getClass().getResource(PATH_ACTION_CREATESELECTION));
//        ACTION_INVERTSELECTION = new ImageIcon(getClass().getResource(PATH_ACTION_INVERTSELECTION));
//        ACTION_INVERTMASK = new ImageIcon(getClass().getResource(PATH_ACTION_INVERTMASK));
//        ACTION_SPECIFYSELECTION = new ImageIcon(getClass().getResource(PATH_ACTION_SPECIFYSELECTION));
//        ACTION_BRUSHSIZE = new ImageIcon(getClass().getResource(PATH_ACTION_BRUSHSIZE));
//        MASK_LOCKED = new ImageIcon(getClass().getResource(PATH_MASK_LOCKED));
//        MASK_UNLOCKED = new ImageIcon(getClass().getResource(PATH_MASK_UNLOCKED));
    	 TOOL_RECTANGLE = new ImageIcon(PATH_TOOL_RECTANGLE);
         TOOL_ROUNDRECTANGLE = new ImageIcon(PATH_TOOL_ROUNDRECTANGLE);
         TOOL_ELLIPSE = new ImageIcon(PATH_TOOL_ELLIPSE);
         TOOL_OVAL = new ImageIcon(PATH_TOOL_OVAL);
         TOOL_POLYGON = new ImageIcon(PATH_TOOL_POLYGON);
         TOOL_FREEHAND = new ImageIcon(PATH_TOOL_FREEHAND);
         TOOL_BRUSH = new ImageIcon(PATH_TOOL_BRUSH);
         ACTION_CREATEMASK = new ImageIcon(PATH_ACTION_CREATEMASK);
         ACTION_CREATESELECTION = new ImageIcon(PATH_ACTION_CREATESELECTION);
         ACTION_INVERTSELECTION = new ImageIcon(PATH_ACTION_INVERTSELECTION);
         ACTION_INVERTMASK = new ImageIcon(PATH_ACTION_INVERTMASK);
         ACTION_SPECIFYSELECTION = new ImageIcon(PATH_ACTION_SPECIFYSELECTION);
         ACTION_BRUSHSIZE = new ImageIcon(PATH_ACTION_BRUSHSIZE);
         MASK_LOCKED = new ImageIcon(PATH_MASK_LOCKED);
         MASK_UNLOCKED = new ImageIcon(PATH_MASK_UNLOCKED);
    }
}
         
