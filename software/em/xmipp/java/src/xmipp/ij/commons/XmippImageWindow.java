package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Label;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;

import javax.swing.SwingUtilities;

import xmipp.ij.commons.XmippMenuBar.IJRequirement;
import xmipp.utils.Params;
import xmipp.utils.XmippWindowUtil;

public class XmippImageWindow extends ImageWindow implements XmippIJWindow
{

	protected XmippMenuBar menu;
	protected ImagePlusLoader ipl;
	protected Label pixelslb;
    private DesignMaskJFrame maskfr;
    protected Params params;
    protected boolean closing = false;
    

	public XmippImageWindow(ImagePlusLoader ipl, Params params)
	{
		this(ipl, ipl.getName(), params);
	}
        
    public XmippImageWindow(ImagePlus imp, Params params)
    {
        this(imp, new XmippImageCanvas(imp), params);
    }
    
    public XmippImageWindow(ImagePlus imp, ImageCanvas canvas, Params params)
	{
		super(imp, canvas);
        this.ipl = new ImagePlusLoader(imp);
        this.params = params;
        XmippApplication.addInstance(true);
        initComponents();
        XmippWindowUtil.setScipionImageIcon(this);
	}


	public XmippImageWindow(ImagePlusLoader ipl, String title, Params params)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.ipl = ipl;
		imp.setTitle(title);
		this.params = params;
		XmippApplication.addInstance(true);
		initComponents();
		XmippWindowUtil.setScipionImageIcon(this);
	}

    protected void initComponents()
    {
    	
        menu = new XmippMenuBar(this);
        if(params != null)//otherwise is picker window
        {
	        setMenuBar(menu);
	        if(XmippApplication.isScipion())
	        {
	            MenuItem createMaskmi = new MenuItem("Design Mask");
	            createMaskmi.addActionListener(new ActionListener() {
	
	                @Override
	                public void actionPerformed(ActionEvent ae) {
	                    loadMaskFrame();
	                }
	            });
	            menu.advancedmn.add(createMaskmi);
	            if(params!= null && params.isMask())
	                loadMaskFrame();
	        }
        }
        pixelslb = new Label("                                                ");
        add(pixelslb);
        

    }

    @Override
	public void openMaskToolbar()
	{
		menu.runCommand("Masks Tool Bar", new IJRequirement[] { IJRequirement.IMAGEJ });
	}

	@Override
	public void loadData()
	{
		getCanvas().loadData(this);
		if(params != null)//otherwise is picker window
			menu.applyFilters();
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
        @Override
	public void windowActivated(WindowEvent e)
	{
		//		if (IJ.isMacintosh())
		//			this.setMenuBar(Menus.getMenuBar());
		
		if (!closed)
		{
			//ic.requestFocus();
			WindowManager.setCurrentWindow(this);
		}
	}

        @Override
	public XmippImageCanvas getCanvas()
	{
		return ((XmippImageCanvas) super.getCanvas());
	}

	public void showPixels(int x, int y, int[] pixels)
	{
        if(pixelslb == null)//avoid error running showPixels before initialization
            return;
        pixelslb.setText(pixelsToString(imp, x, y, pixels));
		
	}
        
    public static String pixelsToString(ImagePlus imp, int x, int y, int[] pixels)
    {
        String text = "";
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
                    text = String.format("x=%s, y=%s, value=%-5s", x, y, value);
                    break;
            case ImagePlus.GRAY32:
                    text = String.format("x=%s, y=%s, value=%.2f", x, y, Float.intBitsToFloat(pixels[0]));
                    break;
            case ImagePlus.COLOR_256:
            case ImagePlus.COLOR_RGB:
                    value =  pixels[0] + "," + pixels[1] + "," + pixels[2];
                    text = String.format("x=%s, y=%s, value=%-15s", x, y, value);
                    break;
        }
        return text;
        
    }
        
        
        @Override
	public void windowClosing(WindowEvent e) {
            
        super.windowClosing(e);
        myClose();
        imp.flush();
	}
        
    public synchronized void myClose()
    {
    	closing = true;
    	XmippApplication.removeInstance(true);
    	if(maskfr != null)
    	{
    		maskfr.setVisible(false);
    		maskfr.dispose();
    		
    	}
    }
    
        
    @Override
    public  boolean close()
    {
        boolean result = super.close();
        myClose();
        return result;
    }
    
   
    public Params getParams()
    {
        return params;
    }


    protected void loadMaskFrame()
    {
        if(maskfr == null || !maskfr.isVisible())
        {
            SwingUtilities.invokeLater(new Runnable() {

                @Override
                public void run() {
                    maskfr = new DesignMaskJFrame(XmippImageWindow.this);
                    XmippWindowUtil.setLocation(0.2f, 0.5f, XmippImageWindow.this);
                    XmippWindowUtil.setLocation(0.6f, 0.5f, maskfr);
                }
            });
            
            imp.getCanvas().addMouseListener(new MouseListener() {

                @Override
                public void mouseClicked(MouseEvent me) {
                }

                @Override
                public void mousePressed(MouseEvent me) {
                }

                @Override
                public void mouseReleased(MouseEvent me) {
                    maskfr.refreshMask();
                }

                @Override
                public void mouseEntered(MouseEvent me) {
                }

                @Override
                public void mouseExited(MouseEvent me) {
                }
            });
            imp.getCanvas().addComponentListener(new java.awt.event.ComponentAdapter()
            {
                @Override
				public void componentResized(ComponentEvent e)
				{
					maskfr.refreshMask();
				}
			});
            
        }
    }

   
	public synchronized boolean isClosing()
	{
		// TODO Auto-generated method stub
		return closing;
	}

	
	
	

}// class XmippImageWindow
