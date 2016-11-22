package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import java.awt.Rectangle;
import java.awt.event.*;
import javax.swing.SwingUtilities;


public class XmippImageCanvas extends ImageCanvas implements MouseWheelListener, KeyListener
{
	
	private boolean invertx;
	private boolean inverty;
    boolean[] keys = new boolean[1024];

	public Tool getTool()
	{

		if (IJ.getInstance() == null || !IJ.getInstance().isVisible())
			return Tool.PICKER;
		return Tool.getTool(IJ.getToolName());
	}

	public XmippImageCanvas(ImagePlus imp)
	{
		super(imp);
		addMouseWheelListener(this);
	

	}

	public void mousePressed(MouseEvent e)
	{
		if (getTool() == Tool.IMAGEJ)
		{
			super.mousePressed(e);
			return;
		}
		if (isDragImage(e))
		{
			customScrollSetup(e);
			return;
		}
		if(e.isControlDown())
		{
			 if (SwingUtilities.isLeftMouseButton(e))
				 zoomIn(e.getX(), e.getY());
			 if (SwingUtilities.isRightMouseButton(e))
				 zoomOut(e.getX(), e.getY());
		}
	}
    public void customScrollSetup(MouseEvent e){

        int x = super.offScreenX(e.getX());
        int y = super.offScreenY(e.getY());

        setupScroll(x, y);

    }

	protected boolean isDragImage(MouseEvent e)
	{
		return (SwingUtilities.isRightMouseButton(e) && !e.isControlDown());
	}

	public void mouseDragged(MouseEvent e)
	{

		if (getTool() == Tool.IMAGEJ)
		{
			super.mouseDragged(e);
			return;
		}
		if (isDragImage(e))
		{
			scroll(e.getX(), e.getY());
			return;
		}
	}

	
       

	public void mouseReleased(MouseEvent e)
	{
		if (getTool() == Tool.IMAGEJ)// do nothing for ImageJ tool is mine is selected
		{
			super.mouseReleased(e);
			return;
		}

	}
	
	public int getXOnImage(int x)
	{
		int x0 = (int) getSrcRect().getX();
		return (int) ((x - x0) * magnification);
	}
	
	public int getYOnImage(int y)
	{
		int y0 = (int) getSrcRect().getY();
		return (int) ((y - y0) * magnification);
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		if (!e.isShiftDown())
			return;

		int x = e.getX();
		int y = e.getY();
		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);


	}

	public void mouseMoved(MouseEvent e)
	{
		int x = offScreenX(e.getX());
		int y = offScreenY(e.getY());
		
		if(getParent() instanceof XmippImageWindow)
		{
			int[] pixels = imp.getPixel(x, y);
            if(pixels != null)
                ((XmippImageWindow)getParent()).showPixels(x, y, pixels);
		}
		if (getTool() == Tool.IMAGEJ)
			super.mouseMoved(e);

		imp.mouseMoved(x, y);
		imp.updateStatusbarValue();
		
	}
	
	public int offScreenX(int x)
	{
		x = super.offScreenX(x);
		if(invertx)
			x = imp.getWidth() - x;
		return x;
	}
	
	public int offScreenY(int y)
	{
		y = super.offScreenY(y);
		if(inverty)
			y = imp.getHeight() - y;
		return y;
	}

	public void loadData(XmippIJWindow xiw)
	{
		double currmagnif = getMagnification();
		Rectangle rect = getSrcRect();
		imp = xiw.getImagePlusLoader().loadImagePlus();
        imp.setTitle(xiw.getImagePlusLoader().getName());
		ImageWindow iw = (ImageWindow) xiw;
		iw.setImage(getImage());
		iw.updateImage(getImage());
                
		setSourceRect(rect);
		double prefmagnif = getPreferredMagnification();
		if (currmagnif < prefmagnif)
			setMagnification(prefmagnif);
		else
			setMagnification(currmagnif);
		setDrawingSize((int) (rect.getWidth() * magnification), (int) (rect.getHeight() * magnification));

		repaint();
		iw.pack();

	}

	public double getPreferredMagnification()
	{
		double magnification = getMagnification();
		int min = 200;
		while (Math.min(getSrcRect().getWidth() * magnification, getSrcRect().getHeight() * magnification) < min)
		{
			magnification = 2 * magnification;
		}
		return magnification;
	}

	public void adjustMagnification()// for micrographs will not happen
	{
		double currmagnif = getMagnification();
		double prefmagnif = getPreferredMagnification();
		if (currmagnif < prefmagnif)
		{
			setMagnification(prefmagnif);
			setDrawingSize((int) (getSrcRect().getWidth() * magnification), (int) (getSrcRect().getHeight() * magnification));
			repaint();
			if (getParent() != null)
			{
				ImageWindow iw = (ImageWindow) getParent();
				iw.pack();
			}
		}
	}
	
	public void setInvertX(boolean value)
	{
		invertx = value;
	}
	
	public void setInvertY(boolean value)
	{
		inverty = value;
	}

    @Override
    public void keyTyped(KeyEvent keyEvent) {

    }

    @Override
    public void keyPressed(KeyEvent keyEvent) {

        keys[keyEvent.getKeyCode()] = true;
    }

    @Override
    public void keyReleased(KeyEvent keyEvent) {
        keys[keyEvent.getKeyCode()] = false;
    }

    protected boolean isKeyPressed(int keyCode){
        if (keyCode < keys.length){
            return keys[keyCode];
        }
        return false;
    }


}
