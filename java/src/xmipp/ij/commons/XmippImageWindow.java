package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageWindow;
import java.awt.Label;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import xmipp.ij.commons.XmippMenuBar.IJRequirement;

public class XmippImageWindow extends ImageWindow implements XmippIJWindow
{

	protected XmippMenuBar menu;
	private ImagePlusLoader ipl;
	private Label pixelslb;

	public XmippImageWindow(ImagePlusLoader ipl)
	{
		this(ipl, ipl.getFileName());
	}

	public XmippImageWindow(ImagePlusLoader ipl, String title)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
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
		pixelslb = new Label("                                                   ");
		add(pixelslb);
		pack();//avoids header bad view
	}

	public void openMaskToolbar()
	{
		menu.runCommand("Masks Tool Bar", new IJRequirement[] { IJRequirement.IMAGEJ });
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
	public void windowActivated(WindowEvent e)
	{
		//		if (IJ.isMacintosh())
		//			this.setMenuBar(Menus.getMenuBar());
		if (IJ.debugMode)
			IJ.write(imp.getTitle() + ": Activated");
		if (!closed)
		{
			//ic.requestFocus();
			WindowManager.setCurrentWindow(this);
		}
	}

	public XmippImageCanvas getCanvas()
	{
		return ((XmippImageCanvas) super.getCanvas());
	}

	public void showPixels(int x, int y, int[] pixels)
	{
		String value = "";
		switch (imp.getType())
		{
		case ImagePlus.GRAY8:
		case ImagePlus.GRAY16:
			double cValue = imp.getCalibration().getCValue(pixels[0]);
			if (cValue == pixels[0])
				value = String.valueOf(pixels[0]);
			else
				value = IJ.d2s(cValue) + " (" + pixels[0] + ")";
			break;
		case ImagePlus.GRAY32:
			value = String.valueOf(Float.intBitsToFloat(pixels[0]));
			break;
		case ImagePlus.COLOR_256:
			value = pixels[0] + "," + pixels[1] + "," + pixels[2];
			break;
		case ImagePlus.COLOR_RGB:
			value = pixels[0] + "," + pixels[1] + "," + pixels[2];
			break;
		
		}
		String text = String.format("x=%s, y=%s, value=%s", x, y, value);
		pixelslb.setText(text);
	}

}// class XmippImageWindow
