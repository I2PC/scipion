package xmipp.ij;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import java.awt.Menu;
import java.awt.MenuBar;
import javax.swing.SwingUtilities;

import xmipp.jni.Filename;


public class XmippImageWindow extends ImageWindow
{
	public XmippImageWindow(ImagePlus imp) {
		this(imp, "");
	}
	
	public XmippImageWindow(ImagePlus imp, String title)
	{
		super(imp, new XmippImageCanvas(imp));
		setTitle(title);
		setMenuBar(new XmippMenuBar());
	}
	
	public static void openImageJ(Tool tool){
		if (IJ.getInstance() == null)
		{
			new ImageJ();
			IJ.run("Install...", "install=" + Filename.getXmippPath("java/src/xmipp/ij/XmippMacros.txt"));
			IJ.setTool(Tool.getTool(tool));
		}		
	}
}//class XmippImageWindow
