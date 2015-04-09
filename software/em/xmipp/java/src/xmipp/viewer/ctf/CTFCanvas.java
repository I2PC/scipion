package xmipp.viewer.ctf;

import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import xmipp.utils.DEBUG;

@SuppressWarnings("serial")
public class CTFCanvas extends ImageCanvas {

	private CTFRecalculateImageWindow master;
	private int width;
	
	public CTFCanvas(ImagePlus imp) {
		super(imp);
	}
	
	public void setMaster(CTFRecalculateImageWindow master){
		this.master = master;
		ImagePlus imp = master.getImagePlus();
		width = imp.getWidth();
	}
	
	public void paint(Graphics g)
	{
		super.paint(g);
		g.setColor(Color.blue);
		paintFreqCircle(g, master.getLowFreq());
		paintFreqCircle(g, master.getHighFreq());
	}
	
	public void paintFreqCircle(Graphics g, double freq){ 
                double size = width * magnification;
		int outer = (int)(size * freq * 2);
		int offset = (int)(size - outer) / 2;
		g.drawOval(offset, offset, outer, outer);
	}

}
