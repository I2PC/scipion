package xmipp.ij;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.io.FileSaver;

import java.io.File;

import xmipp.jni.Filename;


public class XmippImageWindow extends ImageWindow implements XmippIJWindow
{
	
	

	public static void main(String[] args)
	{
		try
		{
			new XmippImageWindow("/home/airen/Coss/Xmipp/BPV_2/InputData/BPV_1386.mrc");
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	private String file;
	
	public XmippImageWindow(String file) throws Exception
	{
		this(XmippImageConverter.loadImage(file), file);
	}
	
	public XmippImageWindow(ImagePlus imp) {
		this(imp, "");
	}
	
	public XmippImageWindow(ImagePlus imp, String title)
	{
		super(imp, new XmippImageCanvas(imp));
		file = imp.getOriginalFileInfo().directory + File.separator + imp.getOriginalFileInfo().fileName;
		setTitle(title);
		setMenuBar(new XmippMenuBar(this));
	}
	
	public static void openImageJ(Tool tool){
		if (IJ.getInstance() == null)
		{
			new ImageJ();
			IJ.run("Install...", "install=" + Filename.getXmippPath("java/src/xmipp/ij/XmippMacros.txt"));
		}		
		boolean recognized = IJ.setTool(Tool.getTool(tool));
		//System.out.println(recognized);
	}

	@Override
	public void loadData()
	{
		
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
		System.out.println(imp.getTitle());
		saveDataAs(imp.getTitle());
		
	}
}//class XmippImageWindow
