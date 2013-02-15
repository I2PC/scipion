package xmipp.ij.commons;

import ij.ImagePlus;
import ij.gui.StackWindow;

import java.awt.Frame;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.WindowEvent;
import xmipp.ij.commons.XmippMenuBar.IJRequirement;

public class XmippStackWindow extends StackWindow implements XmippIJWindow{
	
	private ImagePlusLoader ipl;
	protected XmippMenuBar menu;
	private Window window;

	public XmippStackWindow(Window window, ImagePlusLoader ipl, String title)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.window = window;
		this.ipl = ipl;
		setTitle(title);
		menu = new XmippMenuBar(this);
		setMenuBar(menu);
		((XmippImageCanvas)getCanvas()).adjustMagnification();
	}
	
	public XmippStackWindow(ImagePlusLoader ipl)
	{
		this(null, ipl, ipl.getFileName());
	}
	
	public XmippStackWindow(Window window, ImagePlusLoader ipl)
	{
		this(window, ipl, ipl.getFileName());
	}
	
	public XmippStackWindow(ImagePlusLoader ipl, String title)
	{
		this(null, ipl,title);
	}


	@Override
	public void loadData()
	{
		XmippImageCanvas canvas = (XmippImageCanvas)getCanvas();
		canvas.loadData(this);
		canvas.adjustMagnification();
	}

	public void openMaskToolbar(){
		menu.runCommand("Masks Tool Bar", new IJRequirement[]{IJRequirement.IMAGEJ});
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
		if(window == null)//if I am the main process I can close java
			System.exit(0);
		super.windowClosing(e);
		close();//it works for stack
		if(XmippIJUtil.getXmippImageJ() != null)
			XmippIJUtil.getXmippImageJ().close();
	}

	@Override
	public ImagePlusLoader getImagePlusLoader()
	{
		return ipl;
	}

	@Override
	public boolean isVolume()
	{
		return getImagePlusLoader().isVolume();
	}

	@Override
	public boolean isStack()
	{
		return true;
	}

	


}
