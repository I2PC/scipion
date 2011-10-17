package particlepicker.tiltpair.gui;

import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.swing.SwingUtilities;

import particlepicker.ParticlePickerCanvas;
import particlepicker.ParticlePickerJFrame;
import particlepicker.WindowUtils;
import particlepicker.tiltpair.model.TiltPairPicker;
import particlepicker.tiltpair.model.TiltedParticle;
import particlepicker.tiltpair.model.UntiltedMicrograph;
import particlepicker.tiltpair.model.UntiltedParticle;
import particlepicker.training.model.TrainingParticle;
import xmipp.Particle;

public class UntiltedMicrographCanvas extends ParticlePickerCanvas implements MouseWheelListener
{

	private TiltPairPickerJFrame frame;
	private TiltedMicrographCanvas tiltedcanvas;
	private UntiltedParticle active;
	private TiltPairPicker pppicker;
	private UntiltedMicrograph um;
	private ImageWindow iw;
	private boolean reload = false;
	

	public UntiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getImage());
		this.um = frame.getMicrograph();
	
		this.frame = frame;
		addMouseWheelListener(this);
		this.pppicker = frame.getParticlePairPicker();
		tiltedcanvas = new TiltedMicrographCanvas(frame);
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0, 0, iw);
	}

	public void updateMicrograph()
	{
		this.um = frame.getMicrograph();
		iw.setImage(um.getImage());
		iw.updateImage(um.getImage());
		tiltedcanvas.updateMicrograph();
		active = null;
	}

	public void mouseEntered(MouseEvent e)
	{
		if (frame.getTool() != Tool.PICKER)
		{
			super.mouseEntered(e);
			return;
		}
		setCursor(crosshairCursor);
	}

	public void mouseMoved(MouseEvent e)
	{
		if (frame.getTool() != Tool.PICKER)
		{
			super.mouseMoved(e);
			return;
		}
		setCursor(crosshairCursor);
	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	public void mousePressed(MouseEvent e)
	{
		if (frame.getTool() != Tool.PICKER)
		{
			super.mousePressed(e);
			return;
		}
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (SwingUtilities.isRightMouseButton(e))
		{
			setupScroll(x, y);
			tiltedcanvas.mousePressed(x, y);
		}
		if(active != null && !active.isAdded() && active.getTiltedParticle() != null)
				um.addParticleToAligner(active);
		UntiltedParticle p = um.getParticle(x, y, (int) (frame.getParticleSize()));
		
		if (p != null)
		{
			if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
			{
				if(active!= null && active.equals(p))
					active = um.getParticles().get(um.getParticles().size() - 1);
				um.removeParticle(p);
				frame.updateMicrographsModel();
				if(p.isAdded())
					um.initAligner();
				repaint();
				frame.setChanged(true);
			}
			else if (SwingUtilities.isLeftMouseButton(e))
				setActive(p);
		}
		else if (SwingUtilities.isLeftMouseButton(e) && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			p = new UntiltedParticle(x, y, um, pppicker.getFamily());
			um.addParticle(p);
			setActive(p);
			frame.updateMicrographsModel();
			frame.setChanged(true);
		}
	}
	
	public void setActive(UntiltedParticle up)
	{
		active = up;
		um.setActiveParticle(active);
		TiltedParticle tp = active.getTiltedParticle();
		if(tp == null && um.getAddedCount() >= 4)
			um.setAlignerTiltedParticle(up);
		if(tp != null)
			tiltedcanvas.moveTo(tp);
		repaint();
		tiltedcanvas.repaint();
	}


	/**
	 * Updates particle position and repaints. Sets dragged to null at the end
	 */
	public void mouseReleased(MouseEvent e)
	{
		if (frame.getTool() != Tool.PICKER)
		{
			super.mouseReleased(e);
			return;
		}
		if(reload)
			um.initAligner();
		reload = false;
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	@Override
	public void mouseDragged(MouseEvent e)
	{

		if (frame.getTool() != Tool.PICKER)
		{
			super.mouseDragged(e);
			return;
		}
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (SwingUtilities.isRightMouseButton(e))
		{
			scroll(e.getX(), e.getY());
			tiltedcanvas.mouseDragged(e.getX(), e.getY());
			return;
		}
		if (active != null && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			active.setPosition(x, y);
			if (frame.getParticlesJDialog() != null)
				active.getParticleCanvas(frame).repaint();
			if(active.isAdded())
				reload = true;
		}
		frame.setChanged(true);
		repaint();

	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		tiltedcanvas.setMagnification(magnification);
		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		
		tiltedcanvas.mouseWheelMoved(x, y, rotation);
		
	}

	public TiltedMicrographCanvas getTiltedCanvas()
	{
		return tiltedcanvas;
	}

	public void paint(Graphics g)
	{
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		g2.setColor(frame.getColor());
		int index = 0;

		for (TrainingParticle p : um.getParticles())
		{
			drawShape(g2, p, index == (um.getParticles().size() - 1));
			index++;
		}
		if (um.getActiveParticle() != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, um.getActiveParticle(), true);
		}
	}




	@Override
	public void setActive(TrainingParticle p)
	{
		setActive((UntiltedParticle)p);
	}

	@Override
	public ParticlePickerJFrame getFrame()
	{
		return frame;
	}
	


}
