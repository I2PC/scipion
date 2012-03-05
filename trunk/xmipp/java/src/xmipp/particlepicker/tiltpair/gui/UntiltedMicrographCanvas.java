package xmipp.particlepicker.tiltpair.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePickerCanvas;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.utils.WindowUtil;
import xmipp.utils.XmippMessage;
import xmipp.jni.Particle;

public class UntiltedMicrographCanvas extends ParticlePickerCanvas
{

	private TiltPairPickerJFrame frame;
	private UntiltedParticle active;
	private TiltPairPicker pppicker;
	private UntiltedMicrograph um;
	private ImageWindow iw;
	private boolean reload = false;

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

	@Override
	public Micrograph getMicrograph()
	{
		return um;
	}

	public UntiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getImagePlus());
		this.um = frame.getMicrograph();

		this.frame = frame;
		
		this.pppicker = frame.getParticlePicker();
		iw = new ImageWindow(imp, this);
		WindowUtil.setLocation(0, 0, iw);

	}

	public void updateMicrograph()
	{
		this.um = frame.getMicrograph();
		updateMicrographData();
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
				um.addParticleToAligner(active,true);
			UntiltedParticle p = um.getParticle(x, y, (int) (frame.getParticleSize()));

			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
					removeParticle(p);
				else if (SwingUtilities.isLeftMouseButton(e))
					setActive(p);
			}
			else if (SwingUtilities.isLeftMouseButton(e) && Particle.fits(x, y, frame.getParticleSize(), imp.getWidth(), imp.getHeight()))
				addParticle(x, y);
		}
	}

	public void mouseReleased(MouseEvent e)
	{

		super.mouseReleased(e);
		if (reload)// added particle on matrix has been moved. Matrix changed
					// and tilted particle has to be recalculated
		{
			um.getTiltedMicrograph().removeParticle(active.getTiltedParticle());
			active.setAdded(false);
			um.initAligner();
			um.setAlignerTiltedParticle(active);
			frame.getTiltedCanvas().repaint();
		}
		reload = false;
	}

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
		if (active != null && Particle.fits(x, y, frame.getParticleSize(), imp.getWidth(), imp.getHeight()))
		{
			moveActiveParticle(x, y);
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
		if (frame.drawAngles())
			drawLine(Math.toRadians(um.getUntiltedAngle()), g2);
	}

	private void addParticle(int x, int y)
	{
		try
		{
			Particle tp = um.getAlignerTiltedParticle(x, y);
			ImagePlus tmimage = um.getTiltedMicrograph().getImagePlus();
			if (!Particle.fits(tp.getX(), tp.getY(), pppicker.getFamily().getSize(), tmimage.getWidth(), tmimage.getHeight()))
				throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("Tilted Pair Coordinates"));
			UntiltedParticle p = new UntiltedParticle(x, y, um, pppicker.getFamily());

			um.addParticle(p);

			if (um.getAddedCount() >= 4)
				um.setAlignerTiltedParticle(p);
			setActive(p);
			frame.updateMicrographsModel();
			frame.setChanged(true);
		}
		catch (Exception e)
		{
			JOptionPane.showMessageDialog(this, e.getMessage());
		}

	}

	private void removeParticle(UntiltedParticle p)
	{
		um.removeParticle(p);
		frame.updateMicrographsModel();
		if (active != null && active.equals(p))
		{
			if (!um.getParticles().isEmpty())
				setActive(um.getParticles().get(um.getParticles().size() - 1));
			else
				setActive(null);
		}

		if (p.isAdded())
			um.initAligner();
		repaint();
		frame.getTiltedCanvas().repaint();
		frame.setChanged(true);
	}

	public void setActive(TrainingParticle up)
	{
		active = (UntiltedParticle) up;
		if (active != null)
		{
			TiltedParticle tp = active.getTiltedParticle();
			if (tp != null)
			{
				Rectangle srcrect = frame.getTiltedCanvas().getSrcRect();
				int xrect = (int) ((tp.getX() - srcrect.getX()));
				int yrect = (int) ((tp.getY() - srcrect.getY()));

				if (tp != null && !Particle.fits(xrect, yrect, (int) (tp.getFamily().getSize()), srcrect.width, srcrect.height))
					frame.getTiltedCanvas().moveTo(tp);
			}
		}
		repaint();
		frame.getTiltedCanvas().repaint();
	}

	@Override
	public TrainingParticle getActive()
	{
		return active;
	}

}
