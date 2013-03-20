package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.swing.SwingUtilities;

public class XmippImageCanvas extends ImageCanvas implements MouseWheelListener
{
	protected ImageWindow iw;
	
	public void display()
	{
		if (iw != null && iw.isVisible())
		{
			iw.setImage(getImage());
			iw.updateImage(getImage());
		}
		else
			this.iw = new ImageWindow(getImage(), this);// if you dont provide
														// iw, I init mine
		// iw.maximize();
		iw.pack();
	}


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

		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (isDragImage(e))
		{
			setupScroll(x, y);
			return;
		}

	}

	protected boolean isDragImage(MouseEvent e)
	{
		return SwingUtilities.isRightMouseButton(e) || (SwingUtilities.isLeftMouseButton(e) && e.isControlDown());
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

	public void setZoom(double zoom)
	{
		if (Math.abs(getMagnification() - zoom) <= 0.025)
		{
			if (getMagnification() <= 1.0)
				imp.repaintWindow();
			return;
		}
		if (getMagnification() < zoom)
			zoomIn(0, 0);
		else
			zoomOut(0, 0);

		setZoom(zoom);
	}

	public void mouseReleased(MouseEvent e)
	{
		if (getTool() == Tool.IMAGEJ)// do nothing for ImageJ tool is mine is
										// selected
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
		if (getMagnification() <= 1.0)
			imp.repaintWindow();

	}

	public void mouseMoved(MouseEvent e)
	{

		if (getTool() == Tool.IMAGEJ)
			super.mouseMoved(e);
		int x = offScreenX(e.getX());
		int y = offScreenY(e.getY());
		imp.mouseMoved(x, y);
		imp.updateStatusbarValue();
	}

	public void loadData(XmippIJWindow xiw)
	{
		double currmagnif = getMagnification();
		Rectangle rect = getSrcRect();
		imp = xiw.getImagePlusLoader().loadImagePlus();
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
			System.out.println(magnification);
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



}
