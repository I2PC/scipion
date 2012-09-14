package xmipp.ij.commons;

import ij.ImagePlus;
import ij.gui.StackWindow;

import java.awt.Rectangle;
import java.awt.event.WindowEvent;
import xmipp.ij.commons.XmippMenuBar.IJRequirement;

public class XmippStackWindow extends StackWindow implements XmippIJWindow{
	
	private ImagePlusLoader ipl;
	protected XmippMenuBar menu;

	public XmippStackWindow(ImagePlusLoader ipl, String title)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.ipl = ipl;
		setTitle(title);
		menu = new XmippMenuBar(this);
		setMenuBar(menu);
	}
	
	public XmippStackWindow(ImagePlusLoader ipl)
	{
		this(ipl, ipl.getFileName());
	}

	@Override
	public void loadData()
	{
		XmippImageCanvas canvas = (XmippImageCanvas)getCanvas();
		canvas.loadData(this);
		
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
