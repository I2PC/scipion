package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import xmipp.jni.Particle;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippMessageDialog;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class UntiltedMicrographCanvas extends ParticlePickerCanvas
{

	private TiltPairPickerJFrame frame;
	private UntiltedParticle active;
	private TiltPairPicker pppicker;
	private UntiltedMicrograph um;

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
		super(frame.getMicrograph().getImagePlus(frame.getParticlePicker().getFilters()));
		this.um = frame.getMicrograph();

		this.frame = frame;

		this.pppicker = frame.getParticlePicker();
		um.runImageJFilters(pppicker.getFilters());

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

		if (isDragImage(e))
			frame.getTiltedCanvas().mousePressed(x, y);
		else if (frame.isPickingAvailable(e))
		{
			if (frame.isEraserMode())
			{
				um.removeParticles(x, y);
				active = getLastParticle();
				refresh();

				return;
			}

			if (active != null && !active.isAdded() && active.getTiltedParticle() != null)
				um.addParticleToAligner(active, true);
			UntiltedParticle p = um.getParticle(x, y, (int) (frame.getParticleSize()));

			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
					removeParticle(p);
				else if (SwingUtilities.isLeftMouseButton(e))
					refreshActive(p);
			}
			else if (SwingUtilities.isLeftMouseButton(e))
			{
				if (um.fits(x, y, frame.getParticleSize()))
					addParticle(x, y);
				else
					XmippMessageDialog.showInfo(frame, XmippMessage.getOutOfBoundsMsg(String
							.format("Particle centered at %s, %s with size %s", x, y, frame.getParticlePicker().getSize())));
			}
		}
	}

	protected UntiltedParticle getLastParticle()
	{
		if (um.getParticles().isEmpty())
			return null;
		return um.getParticles().get(um.getParticles().size() - 1);
	}

	@Override
	public void mouseDragged(MouseEvent e)
	{
		super.mouseDragged(e);

		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (isDragImage(e))
		{
			frame.getTiltedCanvas().mouseDragged(e.getX(), e.getY());
			return;
		}
		if (frame.isPickingAvailable(e))
		{
			if (frame.isEraserMode())
			{
				um.removeParticles(x, y);
				active = getLastParticle();
				refresh();

				return;
			}

			if (active != null && um.fits(x, y, frame.getParticleSize()))

			{
				setActiveMoved(true);
				moveActiveParticle(x, y);

			}
		}
		frame.setChanged(true);
		repaint();

	}

	public void mouseReleased(MouseEvent e)
	{

		super.mouseReleased(e);
		int x = e.getX();
		int y = e.getY();
		manageActive(x, y);

	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		super.mouseWheelMoved(e);
		if (!e.isShiftDown())
			return;
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

	@Override
	protected void doCustomPaint(Graphics2D g2)
	{
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
			drawLine(Math.toRadians(um.getUntiltedAngle()), g2);// TODO
																// Auto-generated
																// method stub

	}

	@Override
	public void setMicrograph(Micrograph m)
	{
		um = (UntiltedMicrograph) m;
	}

	private void addParticle(int x, int y)
	{
		try
		{
			Particle tp = um.getAlignerTiltedParticle(x, y);
			if (um.getAddedCount() > UntiltedMicrograph.getAlignmentMin()
					&& !um.getTiltedMicrograph().fits(tp.getX(), tp.getY(), pppicker.getSize()))
				throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("Tilted Pair Coordinates"));
			UntiltedParticle p = new UntiltedParticle(x, y, um, pppicker);

			um.addParticle(p);

			if (um.getAddedCount() >= UntiltedMicrograph.getAlignmentMin())
				um.setAlignerTiltedParticle(p);
			refreshActive(p);
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

		if (active != null && active.equals(p))
		{
			if (!um.getParticles().isEmpty())
				refreshActive(um.getParticles().get(um.getParticles().size() - 1));
			else
				refreshActive(null);
		}

		if (p.isAdded())
			um.initAligner();
		refresh();
		frame.getTiltedCanvas().repaint();
	}

	
	public void refreshActive(Particle up)
	{
		
		if (up != null)
		{
			active = (UntiltedParticle) up;
			TiltedParticle tp = active.getTiltedParticle();
			if (tp != null)
			{
				Rectangle srcrect = frame.getTiltedCanvas().getSrcRect();
				int xrect = (int) ((tp.getX() - srcrect.getX()));
				int yrect = (int) ((tp.getY() - srcrect.getY()));

				if (tp != null && !um.fits(xrect, yrect, pppicker.getSize()))
					frame.getTiltedCanvas().moveTo(tp);
			}
		}
		else
			active = null;
		repaint();
		frame.getTiltedCanvas().repaint();
	}

	@Override
	public TrainingParticle getActive()
	{
		return active;
	}

	protected void manageActive(int x, int y)
	{
		if (!activemoved)
			return;
		if (um.fits(x, y, frame.getParticleSize()))
		{
			moveActiveParticle(x, y);
			um.getTiltedMicrograph().removeParticle(active.getTiltedParticle());
		}
		if (active.isAdded())// added particle on matrix has been moved. Matrix changed and tilted particle has to be recalculated
		{

			active.setAdded(false);
			um.initAligner();			
		}
		um.setAlignerTiltedParticle(active);
		frame.getTiltedCanvas().repaint();
		setActiveMoved(false);
	}



}
