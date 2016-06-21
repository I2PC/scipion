package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;

import javax.swing.SwingUtilities;

import xmipp.jni.Particle;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippMessageDialog;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.model.ManualParticle;

public class UntiltedMicrographCanvas extends ParticlePickerCanvas
{

	private UntiltedParticle active;

	@Override
	public TiltPairPickerJFrame getFrame()
	{
		return (TiltPairPickerJFrame)frame;
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
	public UntiltedMicrograph getMicrograph()
	{
		return (UntiltedMicrograph)micrograph;
	}

	public UntiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame);

	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		if (frame.isPickingAvailable(e))
		{
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());

			if (isDragImage(e))
				getFrame().getTiltedCanvas().mousePressed(x, y);
			else
			{
				if (frame.isEraserMode())
				{
					getMicrograph().removeParticles(x, y);
					active = getLastParticle();
					refresh();
                    getFrame().getTiltedCanvas().repaint();
					return;
				}
				if (active != null && !active.isAdded() && active.getTiltedParticle() != null)
					getMicrograph().addParticleToAligner(active, true);
				UntiltedParticle p = getMicrograph().getParticle(x, y, getParticlePicker().getSize());
				if(p == null && active != null && active.contains(x, y))
					p = active;
				if (p != null)
				{
					if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
						removeParticle(p);
					else if (SwingUtilities.isLeftMouseButton(e))
						refreshActive(p);
				}
				else if (SwingUtilities.isLeftMouseButton(e))
				{
					if (micrograph.fits(x, y, getFrame().getParticleSize()))
						addParticle(x, y);
					else
						XmippMessageDialog.showInfo(frame, XmippMessage.getOutOfBoundsMsg(String
								.format("Particle centered at %s, %s with size %s", x, y, frame.getParticlePicker().getSize())));
				}
			}
		}
	}

	protected UntiltedParticle getLastParticle()
	{
		if (getMicrograph().getParticles().isEmpty())
			return null;
		return getMicrograph().getParticles().get(getMicrograph().getParticles().size() - 1);
	}

	@Override
	public void mouseDragged(MouseEvent e)
	{
		super.mouseDragged(e);

		if (frame.isPickingAvailable(e))
		{
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());
			if (isDragImage(e))
			{
				getFrame().getTiltedCanvas().mouseDragged(e.getX(), e.getY());
				return;
			}

			if (frame.isEraserMode())
			{
				getMicrograph().removeParticles(x, y);
				active = getLastParticle();
				refresh();
				return;
			}

			if (active != null && micrograph.fits(x, y, getFrame().getParticleSize()))

			{
				setActiveMoved(true);
				moveActiveParticle(x, y);

			}
			frame.setChanged(true);
			repaint();
		}

	}

	public void mouseReleased(MouseEvent e)
	{
		super.mouseReleased(e);
		if (frame.isPickingAvailable(e))
		{
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());
			manageActive(x, y);
		}

	}
	
	public void zoomIn(int sx, int sy)
	{
		super.zoomIn(sx, sy);
		
		if(getParticlePicker().getZoom() != getFrame().getTiltedCanvas().getMagnification())
			getFrame().getTiltedCanvas().zoomIn(sx, sy);;
	}
	

	public void zoomOut(int sx, int sy)
	{
		super.zoomOut(sx, sy);
		
		if(getParticlePicker().getZoom() != getFrame().getTiltedCanvas().getMagnification())
			getFrame().getTiltedCanvas().zoomOut(sx, sy);;
	}
        
       

	@Override
	protected void doCustomPaint(Graphics2D g2)
	{
		g2.setColor(getFrame().getColor());
		for (ManualParticle p : getMicrograph().getParticles())
		{
			drawShape(g2, p, false);
		}
		if (active != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, active, true);
		}
		if (getFrame().drawAngles())
			drawLine(Math.toRadians(getMicrograph().getUntiltedAngle()), g2);// TODO
																// Auto-generated
																// method stub

	}

	
	private void addParticle(int x, int y)
	{
		try
		{
			if (active != null && active.getTiltedParticle() == null)
			{
				XmippMessageDialog.showInfo(frame, "Remember to pick tilted particle for each particle");
				return;
			}
			Particle tp = getMicrograph().getAlignerTiltedParticle(x, y);
			if (getMicrograph().getAddedCount() > UntiltedMicrograph.getAlignmentMin() && !getMicrograph().getTiltedMicrograph().fits(tp.getX(), tp.getY(), getParticlePicker().getSize()))
				throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("Tilted particle"));
			UntiltedParticle p = new UntiltedParticle(x, y, getMicrograph(), getParticlePicker());

			getMicrograph().addParticle(p);

			if (tp != null)
				getMicrograph().setAlignerTiltedParticle(p);
			refreshActive(p);
			frame.updateMicrographsModel();
			frame.setChanged(true);
		}
		catch (Exception e)
		{
			XmippDialog.showInfo(frame, e.getMessage());
		}

	}

	public void removeParticle(UntiltedParticle p)
	{
		getMicrograph().removeParticle(p);

		if (active != null && active.equals(p))
		{
			if (!getMicrograph().getParticles().isEmpty())
				refreshActive(getMicrograph().getParticles().get(getMicrograph().getParticles().size() - 1));
			else
				refreshActive(null);
		}

		if (p.isAdded())
			getMicrograph().initAligner();
		refresh();
		getFrame().getTiltedCanvas().repaint();
	}

	public void refreshActive(Particle up)
	{

		if (up != null)
		{
			active = (UntiltedParticle) up;
			TiltedParticle tp = active.getTiltedParticle();
			if (tp != null)
			{
				int x = getFrame().getTiltedCanvas().getXOnImage(tp.getX());
				int y = getFrame().getTiltedCanvas().getYOnImage(tp.getY());

				if (tp != null && !micrograph.fits(x, y, getParticlePicker().getSize()))
					getFrame().getTiltedCanvas().moveTo(tp);
			}
		}
		else
			active = null;
		repaint();
		getFrame().getTiltedCanvas().repaint();
	}

	@Override
	public UntiltedParticle getActive()
	{
		return active;
	}

	protected void manageActive(int x, int y)
	{
		if (!activemoved)
			return;
		boolean hadtilted = active.getTiltedParticle() != null;
		if (micrograph.fits(x, y, getFrame().getParticleSize()))
		{
			moveActiveParticle(x, y);
//			getMicrograph().getTiltedMicrograph().removeParticle(active.getTiltedParticle());
		}
		if (active.isAdded())// added particle on matrix has been moved. Matrix
								// changed and tilted particle has to be
								// recalculated
		{
			active.setAdded(false);
			getMicrograph().initAligner();
		}
//		getMicrograph().setAlignerTiltedParticle(active);
//		if(active.getTiltedParticle() == null && hadtilted)
//			XmippDialog.showInfo(frame, "Tilted particle will be dismissed");
		getFrame().getTiltedCanvas().repaint();
		setActiveMoved(false);
	}

	@Override
	public TiltPairPicker getParticlePicker()
	{
		// TODO Auto-generated method stub
		return (TiltPairPicker)picker;
	}
     


        

}
