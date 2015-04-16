/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.utils;

import java.awt.Image;
import java.io.File;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Dictionary;

import ij.ImagePlus;
import xmipp.jni.Filename;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;

/**
 *
 * @author Juanjo Vega
 */
public class XmippResource {

    public static int DEFAULT_PREVIEW_WIDTH = 80;
    public static int DEFAULT_PREVIEW_HEIGHT = 80;
    
    //Path to resource folder    
    public static String PATH_ICONS;
    //Icon names    
    public final static String MISSING = "missing.png";
    public final static String MISSING2 = "missing2.png";
    public final static String WAIT = "wait.gif";
    public final static String WAIT_MENU = "wait_menu.gif";
    public final static String NORMALIZE = "histogram.png";
    public final static String APPLY_GEO = "rotate.png";
    public final static String ADJUST_COLS = "resize.png";
    public final static String GOTO = "goto.gif";
    public final static String ZOOM = "zoom.png";
    public final static String VIEW_MD = "md_view.gif";
    public final static String VIEW_GALLERY = "table_view.gif";
    
    public static ImageIcon WAIT_ICON;
    public static ImageIcon LOCK_ICON;
    public static ImageIcon DELETE_ICON;
    public static ImageIcon WAIT_MENU_ICON;
    public static ImagePlus MISSING_ITEM;
    public static ImageIcon VIEW_MD_ICON;
    public static ImageIcon VIEW_GALLERY_ICON;
    public static ImageIcon MISSING_ICON;
    
    private static Hashtable<String, ImageIcon> iconsPool;
    

    static { //this block will be called just once, at class initialization
    	try {
    		PATH_ICONS = Filename.getXmippPath();//"/home/josem/xmipp_current/resources/";
    	}
    	catch (Exception e){
    		PATH_ICONS = ".";
    	}
    	DEBUG.printMessage(PATH_ICONS);
    	iconsPool = new Hashtable<String, ImageIcon>();
    	
    	WAIT_ICON = getIcon(WAIT);
    	LOCK_ICON = getIcon("locked.png");
    	DELETE_ICON = getIcon("delete_item.png");
    	WAIT_MENU_ICON = getIcon(WAIT_MENU);
    	VIEW_MD_ICON = getIcon(VIEW_MD);
    	VIEW_GALLERY_ICON = getIcon(VIEW_GALLERY);
    	
    	// For missing items.
    	MISSING_ICON = getIcon(MISSING2);
    	MISSING_ITEM = new ImagePlus("X", MISSING_ICON.getImage());
    }
    
    // Create an icon using the xmipp resource path
    public static ImageIcon getIcon(String name){
    	ImageIcon ii;
    	if (iconsPool.containsKey(name))
    		ii = iconsPool.get(name);
    	else {
    		ii = new ImageIcon(String.format("%s%sresources%s%s", PATH_ICONS, File.separator, File.separator, name));
    		iconsPool.put(name, ii);
    	}
    	return ii;
    }
    
    public static String getPath(String resource)
    {
    	try
		{
			return String.format("%s%sresources%s%s", Filename.getXmippPath(), File.separator, File.separator, resource);
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}
    }

	public static Image getImage(String name)
	{
		try
		{
			return ImageIO.read(new File(getPath(name)));
		}
		catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}
	}
}
