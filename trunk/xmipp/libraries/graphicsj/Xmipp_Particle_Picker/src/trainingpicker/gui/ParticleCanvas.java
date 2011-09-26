package trainingpicker.gui;

import ij.ImagePlus;
import ij.gui.ImageCanvas;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;

public class ParticleCanvas extends ImageCanvas {

	public ParticleCanvas(ImagePlus imp) {
		super(imp);
	}

	@Override
	public void mousePressed(MouseEvent e) {
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		
		System.out.printf("%s %s\n", x, y);
	}
	
	public void paint(Graphics g)
	{
		super.paint(g);
		Graphics2D g2 = (Graphics2D)g;
		g2.drawRect(0, 0, imp.getWidth(), imp.getHeight());
	}

	
	


}
