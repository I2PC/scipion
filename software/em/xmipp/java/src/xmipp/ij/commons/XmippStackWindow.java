package xmipp.ij.commons;

import ij.gui.StackWindow;

import java.awt.Window;
import java.awt.event.WindowEvent;

import xmipp.ij.commons.XmippMenuBar.IJRequirement;
import xmipp.utils.Params;

public class XmippStackWindow extends StackWindow implements XmippIJWindow{
	
	private ImagePlusLoader ipl;
	protected XmippMenuBar menu;
	protected Params params;

	public XmippStackWindow(ImagePlusLoader ipl, String title, Params params)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.ipl = ipl;
		imp.setTitle(title);
		menu = new XmippMenuBar(this);
		setMenuBar(menu);
		//getCanvas().adjustMagnification();
		XmippApplication.addInstance(true);
		this.params = params;
		
	}
	
	public XmippStackWindow(ImagePlusLoader ipl, Params params)
	{
		this(ipl, ipl.getName(), params);
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

	public Params getParams()
    {
        return params;
    }


}
