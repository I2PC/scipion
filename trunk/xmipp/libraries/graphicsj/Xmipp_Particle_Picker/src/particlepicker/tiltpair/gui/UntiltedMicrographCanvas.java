package particlepicker.tiltpair.gui;

import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
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

public class UntiltedMicrographCanvas extends ParticlePickerCanvas
{

	private TiltPairPickerJFrame frame;
	private UntiltedParticle active;
	private TiltPairPicker pppicker;
	private UntiltedMicrograph um;
	private ImageWindow iw;
	private boolean reload = false;
	private boolean drawalpha;

	public UntiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getImagePlus());
		this.um = frame.getMicrograph();

		this.frame = frame;
		addMouseWheelListener(this);
		this.pppicker = frame.getParticlePicker();
		drawalpha = false;
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0, 0, iw);
		
	}

	public void updateMicrograph()
	{
		this.um = frame.getMicrograph();
		iw.setImage(um.getImagePlus());
		iw.updateImage(um.getImagePlus());
		if (!um.getParticles().isEmpty())
			setActive(um.getParticles().get(um.getParticles().size() - 1));
		else
			setActive(null);
	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (SwingUtilities.isRightMouseButton(e))
			frame.getTiltedCanvas().mousePressed(x, y);
		
		if (frame.isPickingAvailable(e))
		{
			if (active != null && !active.isAdded() && active.getTiltedParticle() != null)
				um.addParticleToAligner(active);
			UntiltedParticle p = um.getParticle(x, y, (int) (frame.getParticleSize()));

			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
				{
					if (active != null && active.equals(p))
						setActive(um.getParticles().get(um.getParticles().size() - 1));
					um.removeParticle(p);
					frame.updateMicrographsModel();
					if (p.isAdded())
						um.initAligner();
					repaint();
					frame.getTiltedCanvas().repaint();
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
	}

	public void setActive(UntiltedParticle up)
	{
		active = up;
		if (active != null)
		{
			TiltedParticle tp = active.getTiltedParticle();
			if (tp == null && um.getAddedCount() >= 4)
				um.setAlignerTiltedParticle(up);
			if (tp != null)
				frame.getTiltedCanvas().moveTo(tp);
		}
		repaint();
		frame.getTiltedCanvas().repaint();
	}

	/**
	 * Updates particle position and repaints. Sets dragged to null at the end
	 */
	public void mouseReleased(MouseEvent e)
	{
		super.mouseReleased(e);
		if (reload)
			um.initAligner();
		reload = false;
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	@Override
	public void mouseDragged(MouseEvent e)
	{
		super.mouseDragged(e);

		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (SwingUtilities.isRightMouseButton(e))
		{
			frame.getTiltedCanvas().mouseDragged(e.getX(), e.getY());
			return;
		}
		if (active != null && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			active.setPosition(x, y);
			if (frame.getParticlesJDialog() != null)
				active.getParticleCanvas(frame).repaint();
			if (active.isAdded())
				reload = true;
		}
		frame.setChanged(true);
		repaint();

	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int x = e.getX();
		int y = e.getY();
		frame.getTiltedCanvas().setMagnification(magnification);
		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		if (getMagnification() <= 1.0)
			imp.repaintWindow();
		frame.getTiltedCanvas().mouseWheelMoved(x, y, rotation);
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
		if (active != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, active, true);
		}
		drawalpha = true;
		if(drawalpha)
		{
			int [] alphas = um.getAlphas();
			double alpha = Math.toRadians(alphas[0]);
			drawLine(alpha, g2);
		}
	}
	
	
	

	@Override
	public void setActive(TrainingParticle p)
	{
		setActive((UntiltedParticle) p);
	}

	@Override
	public ParticlePickerJFrame getFrame()
	{
		return frame;
	}

	public TiltedParticle getActiveTiltedParticle()
	{
		if (active == null)
			return null;
		return active.getTiltedParticle();
	}

	public UntiltedParticle getActiveParticle()
	{
		return active;
	}

	public boolean hasActiveParticle()
	{
		return active != null;
	}

}
