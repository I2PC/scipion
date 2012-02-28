package xmipp.ij;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Menu;
import java.awt.MenuBar;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import xmipp.particlepicker.Tool;
import xmipp.utils.WindowUtil;

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
					ImagePlus imp = new ImagePlus("/home/airen/CellClassifier_/images/8254.tif");
					XmippImageWindow frame = new XmippImageWindow(new XmippImageCanvas(imp));
					
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

	private MenuBar mb;

	public XmippImageWindow(ImageCanvas canvas)
	{
		super(canvas.getImage(), canvas);
		
		
		initComponents();
	}
	
	private void initComponents()
	{
		setTitle("Xmipp Micrograph Viewer");
		initMenuBar();
		setMenuBar(mb);
		
	}
	
	private void initMenuBar()
	{
		mb = new MenuBar();
		mb.add(new Menu("File"));
	}

	
	
}
