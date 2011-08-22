package gui;

import ij.ImagePlus;
import ij.gui.ImageCanvas;

import java.awt.event.MouseEvent;

public class ParticleImageCanvas extends ImageCanvas {

	public ParticleImageCanvas(ImagePlus imp) {
		super(imp);
		
		// TODO Auto-generated constructor stub
	}

	@Override
	public void mousePressed(MouseEvent e) {
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		System.out.printf("%s %s\n", x, y);
	}

	
	


}
