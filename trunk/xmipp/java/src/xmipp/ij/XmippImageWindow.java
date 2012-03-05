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
	
	
	public static void main(String[] args)
	{
		try
		{
			ImagePlus imp = XmippImageConverter.loadImage("/home/airen/Coss/Xmipp/BPV_2/InputData/BPV_1386.mrc");
			new XmippImageWindow(imp);
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
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
		}		
		boolean recognized = IJ.setTool(Tool.getTool(tool));
		System.out.println(recognized);
	}
}//class XmippImageWindow
