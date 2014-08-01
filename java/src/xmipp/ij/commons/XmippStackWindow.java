package xmipp.ij.commons;

import ij.gui.StackWindow;

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
		imp.setTitle(title);
		menu = new XmippMenuBar(this);
		setMenuBar(menu);
		//getCanvas().adjustMagnification();
		XmippApplication.addInstance(true);
		
	}
	
	public XmippStackWindow(ImagePlusLoader ipl)
	{
		this(null, ipl, ipl.getName());
	}
	
	public XmippStackWindow(Window window, ImagePlusLoader ipl)
	{
		this(window, ipl, ipl.getName());
	}
	
	public XmippStackWindow(ImagePlusLoader ipl, String title)
	{
		this(null, ipl,title);
	}


	@Override
	public void loadData()
	{
		getCanvas().loadData(this);
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
	
	public XmippImageCanvas getCanvas()
	{
		return ((XmippImageCanvas)super.getCanvas());
	}
	
	@Override
	public void windowClosing(WindowEvent e) {
 
		 super.windowClosing(e);
                 XmippApplication.removeInstance(true);
               
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
