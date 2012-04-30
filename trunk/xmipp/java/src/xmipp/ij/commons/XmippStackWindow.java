package xmipp.ij.commons;

import ij.ImagePlus;
import ij.gui.StackWindow;

import java.awt.event.WindowEvent;

public class XmippStackWindow extends StackWindow implements XmippIJWindow{
	
	private ImagePlusLoader ipl;


	public XmippStackWindow(ImagePlusLoader ipl, String title)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.ipl = ipl;
		setTitle(title);
		setMenuBar(new XmippMenuBar(this));
	}
	
	public XmippStackWindow(ImagePlusLoader ipl)
	{
		this(ipl, ipl.getFileName());
	}

	@Override
	public void loadData()
	{
		ImagePlus imp = ipl.loadImagePlus();
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
		close();//it works for stack
		if(XmippIJUtil.getXmippImageJ() != null)
			XmippIJUtil.getXmippImageJ().close();
	}

	@Override
	public ImagePlusLoader getImagePlusLoader()
	{
		return ipl;
	}


}
