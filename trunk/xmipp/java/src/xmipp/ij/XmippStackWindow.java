package xmipp.ij;

import xmipp.jni.Filename;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.StackWindow;

public class XmippStackWindow extends StackWindow implements XmippIJWindow{
	public XmippStackWindow(ImagePlus imp) {
		this(imp, "");
	}
	
	public XmippStackWindow(ImagePlus imp, String title)
	{
		super(imp, new XmippImageCanvas(imp));
		setTitle(title);
		setMenuBar(new XmippMenuBar(this));
	}
	
	public static void openImageJ(Tool tool){
		if (IJ.getInstance() == null)
		{
			new ImageJ();
			IJ.run("Install...", "install=" + Filename.getXmippPath("java/src/xmipp/ij/XmippMacros.txt"));
			IJ.setTool(Tool.getTool(tool));
		}		
	}

	@Override
	public void loadData()
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void saveData()
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void saveDataAs(String file)
	{
		// TODO Auto-generated method stub
		
	}
}
