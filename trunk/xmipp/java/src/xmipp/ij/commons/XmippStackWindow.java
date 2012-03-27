package xmipp.ij.commons;

import java.awt.event.WindowEvent;
import java.io.File;

import xmipp.jni.Filename;
import xmipp.utils.DEBUG;
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

	@Override
	public void loadData()
	{
		String file = imp.getOriginalFileInfo().directory + File.separator + imp.getOriginalFileInfo().fileName;
		ImagePlus imp = new ImagePlus(file);
		setImage(imp);//second alone does not work
		updateImage(imp);//first one alone does not work
	}

	

	@Override
	public void saveDataAs(String file) throws Exception
	{
		XmippImageConverter.writeImagePlus(imp, file);
	}

	@Override
	public void saveData() throws Exception
	{
		saveDataAs(imp.getTitle());
		
	}
	
	@Override
	public void windowClosing(WindowEvent e) {
		super.windowClosing(e);
		if(XmippIJUtil.getXmippImageJ() != null)
			XmippIJUtil.getXmippImageJ().close();
		DEBUG.printMessage("kkk");
	}

	@Override
	public String getImageFilePath()
	{
		return imp.getOriginalFileInfo().directory + File.separator + imp.getOriginalFileInfo().fileName;
	}
}
