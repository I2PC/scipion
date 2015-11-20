package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.util.List;

import javax.swing.SwingUtilities;

import xmipp.jni.Particle;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippMessageDialog;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.model.ManualParticle;

public class TiltedMicrographCanvas extends ParticlePickerCanvas
{

	private UntiltedMicrograph um;
	private UntiltedMicrographCanvas uc;

	public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame, frame.getMicrograph().getTiltedMicrograph());
		this.um = frame.getMicrograph();
		this.uc = (UntiltedMicrographCanvas) frame.getCanvas();
		//XmippWindowUtil.setLocation(0.7f, 0, iw);
	}
	
	

	public void updateMicrograph()
	{
		um = getFrame().getMicrograph();
		TiltedMicrograph m = getFrame().getMicrograph().getTiltedMicrograph();
		updateMicrograph(m);
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

			if (frame.isPickingAvailable(e))
			{
				TiltedParticle p = um.getTiltedMicrograph().getParticle(x, y, (int) (getFrame().getParticleSize()));
				if (p != null)
				{
					if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
						uc.removeParticle(p.getUntiltedParticle());
					else if (SwingUtilities.isLeftMouseButton(e))
						refreshActive(p);
				}
				else if (SwingUtilities.isLeftMouseButton(e) && um.fits(x, y, getFrame().getParticleSize()))
				{
					TiltedParticle tp = getMicrograph().getParticle(x, y, getParticlePicker().getSize());
					//associate tilted particle for first cases
					if (um.getAddedCount() < UntiltedMicrograph.getAlignmentMin())
					{
						UntiltedParticle uactive = uc.getActiveParticle();
						if (uactive != null && uactive.getTiltedParticle() == null)
						{
							p = new TiltedParticle(x, y, uactive);
							uactive.setTiltedParticle(p);
							um.getTiltedMicrograph().addParticle(p);
							um.addParticleToAligner(uactive, true);
							frame.updateMicrographsModel();
						}
					}
					else if(tp == null)//generate untilted particle
						addParticle(x, y);
				}
				frame.setChanged(true);
				repaint();
			}
		}
	}

	public void mousePressed(int x, int y)
	{
		setupScroll(x, y);
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	public void mouseDragged(int x, int y)
	{
		scroll(x, y);
	}

	@Override
	public void mouseDragged(MouseEvent e)
	{

		super.mouseDragged(e);
		if (frame.isPickingAvailable(e))
		{
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());

			if (uc.getActive().getTiltedParticle() != null && um.fits(x, y, getFrame().getParticleSize()))
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
			super.mouseReleased(e);
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());
			manageActive(x, y);
		}
	}
        
	public void zoomIn(int sx, int sy)
	{
		super.zoomIn(sx, sy);
		
		if(getParticlePicker().getZoom() != uc.getMagnification())
			uc.zoomIn(sx, sy);	
	}
	
	public void zoomOut(int sx, int sy)
	{
		super.zoomOut(sx, sy);
		
		if(getParticlePicker().getZoom() != uc.getMagnification())
			uc.zoomOut(sx, sy);	
	}

	@Override
	public void refreshActive(Particle p)
	{
		if (p != null)
			frame.getCanvas().refreshActive(((TiltedParticle) p).getUntiltedParticle());
	}

	@Override
	public TiltPairPickerJFrame getFrame()
	{
		return (TiltPairPickerJFrame)frame;
	}

	@Override
	public TiltedMicrograph getMicrograph()
	{
		return (TiltedMicrograph)micrograph;
	}

	@Override
	public TiltedParticle getActive()
	{
		return uc.getActiveTiltedParticle();
	}

	@Override
	protected void doCustomPaint(Graphics2D g2)
	{
		g2.setColor(getFrame().getColor());
		List<TiltedParticle> particles = um.getTiltedMicrograph().getParticles();
		for (TiltedParticle p : particles)
		{
			drawShape(g2, p, false);
		}

		if (uc.getActiveTiltedParticle() != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, uc.getActiveParticle().getTiltedParticle(), true);
		}
		if (getFrame().drawAngles())
			drawLine(Math.toRadians(um.getTiltedAngle()), g2);

	}

	@Override
	protected TiltedParticle getLastParticle()
	{
		return uc.getLastParticle().getTiltedParticle();
	}

	protected void manageActive(int x, int y)
	{
		if (!activemoved)
			return;
		if (uc.getActive().isAdded())
			um.initAligner();
		setActiveMoved(false);
	}

	

	@Override
	public TiltPairPicker getParticlePicker()
	{
		// TODO Auto-generated method stub
		return (TiltPairPicker)picker;
	}
	
	private void addParticle(int x, int y)
	{
		try
		{
			
			Particle p = um.getAlignerUntiltedParticle(x, y);
			if (!um.fits(p.getX(), p.getY(), getParticlePicker().getSize()))
				throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("Untilted particle"));
			UntiltedParticle up = new UntiltedParticle(p.getX(), p.getY(), um, getParticlePicker());
			TiltedParticle tp = new TiltedParticle(x, y, up);
			up.setTiltedParticle(tp);
			um.addParticle(up);
			getMicrograph().addParticle(tp);
			refreshActive(tp);
			frame.updateMicrographsModel();
			frame.setChanged(true);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			XmippDialog.showInfo(frame, e.getMessage());
		}

	}


}
