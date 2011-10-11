package particlepicker.training.gui;

import ij.gui.ImageCanvas;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import particlepicker.ParticlePickerJFrame;
import particlepicker.tiltpair.gui.TiltPairPickerJFrame;
import particlepicker.tiltpair.model.TiltedMicrograph;
import particlepicker.tiltpair.model.TiltedParticle;
import particlepicker.tiltpair.model.UntiltedMicrograph;
import particlepicker.tiltpair.model.UntiltedParticle;
import particlepicker.training.model.TrainingMicrograph;
import particlepicker.training.model.TrainingParticle;

public class ParticleCanvas extends ImageCanvas implements MouseMotionListener, MouseListener
{

	private TrainingParticle particle;
	private int size;
	private int lastx, lasty;
	private boolean dragged;
	private ParticlePickerJFrame frame;
	private int side;

	public ParticleCanvas(TrainingParticle particle, ParticlePickerJFrame frame)
	{
		super(particle.getMicrograph().getImage());
		this.particle = particle;
		this.frame = frame;
		setMagnification(frame.getMagnification());
		this.size = (int)(frame.getFamily().getSize());
		side = (int)(size * magnification);
		setDrawingSize(side, side);
		addMouseMotionListener(this);
		addMouseListener(this);
	}
	
	

	public void paint(Graphics g)
	{
		Rectangle source = new Rectangle(particle.getX0(), particle.getY0(), size, size); 
		setSourceRect(source);
		super.paint(g);
		g.setColor(Color.black);
		g.drawRect(0, 0, side - 1, side - 1);
	}

	@Override
	public void mouseDragged(MouseEvent e)
	{
		if (dragged)
		{
			int movex = lastx - e.getX();
			int movey = lasty - e.getY();
			int x = particle.getX() + movex;
			int y = particle.getY() + movey;
			try
			{
				particle.setPosition(x, y);
				repaint();
				frame.getCanvas().repaint();
			}
			catch (Exception ex)
			{
				JOptionPane.showMessageDialog(this, ex.getMessage());
			}
		}
		lastx = e.getX();
		lasty = e.getY();

	}
	

	@Override
	public void mousePressed(MouseEvent e)
	{
		dragged = true;
		lastx = e.getX();
		lasty = e.getY();
		
		if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown()) 
		{
			if(frame instanceof TrainingPickerJFrame)
				((TrainingMicrograph)frame.getMicrograph()).removeParticle(particle, ((TrainingPickerJFrame)frame).getParticlePicker());
			else if (frame instanceof TiltPairPickerJFrame)
			{
				if(frame.getMicrograph() instanceof UntiltedMicrograph)
					((UntiltedMicrograph)frame.getMicrograph()).removeParticle((UntiltedParticle)particle);
				else if(frame.getMicrograph() instanceof TiltedMicrograph)
					((TiltedMicrograph)frame.getMicrograph()).removeParticle((TiltedParticle)particle);
			}
			frame.getCanvas().repaint();
			frame.updateMicrographsModel();
			frame.loadParticles();
		}
		else
			frame.getCanvas().moveTo(particle);
			
	}

	@Override
	public void mouseReleased(MouseEvent e)
	{
		dragged = false;

	}

}
