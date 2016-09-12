package xmipp.viewer.particlepicker;

import ij.gui.ImageCanvas;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.*;
import javax.swing.SwingUtilities;
import xmipp.utils.XmippDialog;
import xmipp.viewer.particlepicker.tiltpair.gui.TiltPairPickerJFrame;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerJFrame;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

public class ParticleCanvas extends ImageCanvas implements MouseMotionListener, MouseListener
{

	private PickerParticle particle;
	private int size;
	private int lastx, lasty;
	private boolean dragged;
	private ParticlePickerJFrame frame;
	private ParticlePickerCanvas canvas;
	private int side;


    public ParticleCanvas(PickerParticle pickerParticle, ParticlePickerJFrame frame)
	{
		super(pickerParticle.getMicrograph().getImagePlus());
		this.particle = pickerParticle;
		this.frame = frame;
		this.canvas = (pickerParticle instanceof TiltedParticle) ? ((TiltPairPickerJFrame) frame).getTiltedCanvas() : frame.getCanvas();

		this.size = (int) (frame.getParticlePicker().getSize());
		side = frame.getSide();
		setMagnification((float) side / size);
		setDrawingSize(side, side);
		addMouseMotionListener(this);
		addMouseListener(this);
	}

	public void paint(Graphics g)
	{
		Rectangle source = new Rectangle(particle.getX0(size), particle.getY0(size), size, size);
		setSourceRect(source);
		super.paint(g);
		g.setColor(Color.white);
		g.drawRect(0, 0, side - 1, side - 1);
                if(particle instanceof UntiltedParticle)
                {
                    String type = "U";
                    g.drawString(type, 5, side - 2);
                }
                else if(particle instanceof TiltedParticle)
                {
                    String type = "T";
                    g.drawString(type, 5, side - 2);
                }
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
				canvas.repaint();
				frame.setChanged(true);
			}
			catch (Exception ex)
			{
				canvas.repaint();
				repaint();
				XmippDialog.showInfo(frame, ex.getMessage());
				
			}
		}
		lastx = e.getX();
		lasty = e.getY();
	}

	@Override
	public void mousePressed(MouseEvent e)
	{
		if (frame.isPickingAvailable(e))
		{
			dragged = true;
			lastx = e.getX();
			lasty = e.getY();

			if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown()) {
                if (frame instanceof SupervisedPickerJFrame)
                    ((SupervisedPickerMicrograph) frame.getMicrograph()).removeParticle(particle, ((SupervisedPickerJFrame) frame).getParticlePicker());
                else if (frame instanceof TiltPairPickerJFrame) {
                    if (frame.getMicrograph() instanceof UntiltedMicrograph)
                        ((UntiltedMicrograph) frame.getMicrograph()).removeParticle((UntiltedParticle) particle);
                    else if (frame.getMicrograph() instanceof TiltedMicrograph)
                        ((TiltedMicrograph) frame.getMicrograph()).removeParticle((TiltedParticle) particle);
                }
                canvas.repaint();
                frame.updateMicrographsModel();
                frame.loadParticles(false);
            }else
			{
				canvas.moveTo(particle);
				canvas.refreshActive(particle);
			}
		}

	}

	@Override
	public void mouseReleased(MouseEvent e)
	{
		dragged = false;

	}
}
