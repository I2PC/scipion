package xmipp.ij;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import java.awt.Menu;
import java.awt.MenuBar;
import javax.swing.SwingUtilities;

import xmipp.jni.Filename;


public class XmippImageWindow extends ImageWindow
{
	
	public static void main(String[] args)
	{
		SwingUtilities.invokeLater(new Runnable()
		{

			@Override
			public void run()
			{
				try
				{
					ImagePlus imp = XmippImageConverter.loadImage("/home/airen/Coss/Xmipp/BPV_2/InputData/BPV_1387.mrc");
					XmippImageWindow frame = new XmippImageWindow(imp);
					
				}
				catch (Exception e)
				{
					// TODO Auto-generated catch block
					e.printStackTrace();
					System.exit(0);
				}
			}
		});
	}

	private XmippMenuBar mb;

	public XmippImageWindow(ImagePlus imp) {
		this(imp, "");
	}
	
	public XmippImageWindow(ImagePlus imp, String title)
	{
		super(imp, new XmippImageCanvas(imp));
		setTitle(title);
		initComponents();
		this.setVisible(false); //doesn't show by default
	}
	
	private void initComponents()
	{
		mb = new XmippMenuBar();
		setMenuBar(mb);
	}
	
	public static void openImageJ(Tool tool){
		if (IJ.getInstance() == null)
		{
			new ImageJ();
			IJ.run("Install...", "install=" + Filename.getXmippPath("java/src/xmipp/ij/XmippMacros.txt"));
			IJ.setTool(Tool.getTool(tool));
		}		
	}

	
	
}
