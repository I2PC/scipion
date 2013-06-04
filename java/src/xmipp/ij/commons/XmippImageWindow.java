package xmipp.ij.commons;

import ij.IJ;
import ij.WindowManager;
import ij.gui.ImageWindow;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Window;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import xmipp.ij.commons.XmippMenuBar.IJRequirement;
import xmipp.utils.DEBUG;


public class XmippImageWindow extends ImageWindow implements XmippIJWindow
{

	protected XmippMenuBar menu;
	private ImagePlusLoader ipl;
	private Window window;


	public XmippImageWindow(Window window, ImagePlusLoader ipl)
	{
		this(window, ipl, ipl.getFileName());
	}


	public XmippImageWindow(Window window, ImagePlusLoader ipl, String title)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.window = window;
		this.ipl = ipl;
		setTitle(title);
		menu = new XmippMenuBar(this);
		setMenuBar(menu);		
		XmippApplication.addInstance();
		addWindowListener(new WindowAdapter()
		{
			@Override
			public void windowClosing(WindowEvent arg0)
			{
				XmippApplication.removeInstance();
			}
		});


	}
	
	public void openMaskToolbar(){
		menu.runCommand("Masks Tool Bar", new IJRequirement[]{IJRequirement.IMAGEJ});
	}

	
	@Override
	public void loadData()
	{
		getCanvas().loadData(this);
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
	

	
	public ImagePlusLoader getImagePlusLoader()
	{
		return ipl;
	}


	@Override
	public boolean isVolume()
	{
		return false;
	}


	@Override
	public boolean isStack()
	{
		return false;
	}
	
	//overwriting ImageJ event to avoid switching menu
	public void windowActivated(WindowEvent e) {
//		if (IJ.isMacintosh())
//			this.setMenuBar(Menus.getMenuBar());
		if (IJ.debugMode) IJ.write(imp.getTitle() + ": Activated");
		if (!closed) {
			//ic.requestFocus();
			WindowManager.setCurrentWindow(this);
		}
	}
	
	public XmippImageCanvas getCanvas()
	{
		return ((XmippImageCanvas)super.getCanvas());
	}



	

}// class XmippImageWindow
