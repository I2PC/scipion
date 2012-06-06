package xmipp.ij.commons;

import ij.ImagePlus;
import ij.gui.ImageWindow;
import java.awt.event.WindowEvent;

public class XmippImageWindow extends ImageWindow implements XmippIJWindow
{

	

	public static void main(String[] args)
	{
		try
		{
			// openImageJ(Tool.VIEWER);
			//XmippStackWindow w = new XmippStackWindow(new ImagePlusLoader("/home/airen/hand.vol"));
			XmippImageWindow w = new XmippImageWindow(new ImagePlusLoader("/home/airen/xprojects/1/Micrographs/KLH_Dataset_I_Training_0001.mrc"));
			// IJ.open( "/home/airen/Coss/Xmipp/BPV_2/InputData/BPV_1386.mrc");

		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	private ImagePlusLoader ipl;


	public XmippImageWindow(ImagePlusLoader ipl)
	{
		this(ipl, ipl.getFileName());
		
	}


	public XmippImageWindow(ImagePlusLoader ipl, String title)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.ipl = ipl;
		setTitle(title);
		setMenuBar(new XmippMenuBar(this));
		
	}

	
	@Override
	public void loadData()
	{
		try
		{
				ImagePlus imp = ipl.loadImagePlus();
				setImage(imp);// second alone does not work
				updateImage(imp);// first one alone does not work
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
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
	}


	
	public ImagePlusLoader getImagePlusLoader()
	{
		return ipl;
	}


	@Override
	public boolean isVolume()
	{
		// TODO Auto-generated method stub
		return false;
	}


	@Override
	public boolean isStack()
	{
		
		return false;
	}
}// class XmippImageWindow
