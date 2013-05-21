package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.util.List;
import javax.swing.SwingUtilities;
import xmipp.jni.Particle;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class TiltedMicrographCanvas extends ParticlePickerCanvas
{

	private TiltPairPickerJFrame frame;
	private UntiltedMicrograph um;
	private UntiltedMicrographCanvas uc;
	private TiltedParticle active;


	public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getTiltedMicrograph().getImagePlus(frame.getParticlePicker().getFilters()));
		this.um = frame.getMicrograph();
		this.frame = frame;
		this.uc = (UntiltedMicrographCanvas)frame.getCanvas();
		//XmippWindowUtil.setLocation(0.7f, 0, iw);
		um.getTiltedMicrograph().runImageJFilters(frame.getParticlePicker().getFilters());
	}

	public void updateMicrograph()
	{
		UntiltedMicrograph m = frame.getMicrograph();
		setMicrograph(m);
		imp = m.getTiltedMicrograph().getImagePlus(getFrame().getParticlePicker().getFilters());
		m.getTiltedMicrograph().runImageJFilters(getFrame().getParticlePicker().getFilters());
		refreshActive(null);
		
	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	

	public void mouseWheelMoved(int x, int y, int rotation)
	{
		getIw().setSize(uc.getIw().getSize());
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		if (getMagnification() <= 1.0)
			imp.repaintWindow();

	}


	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (frame.isPickingAvailable(e))
		{
			TiltedParticle p = um.getTiltedMicrograph().getParticle(x, y, (int) (frame.getParticleSize()));
			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
				{
					um.removeParticle(p.getUntiltedParticle());
					frame.updateMicrographsModel();
					frame.getCanvas().repaint();
				}
				else if (SwingUtilities.isLeftMouseButton(e))
					active = p;
			}
			else if (uc.hasActiveParticle() && SwingUtilities.isLeftMouseButton(e) && um.fits(x, y, frame.getParticleSize()))
			{
				UntiltedParticle uactive = uc.getActiveParticle();
				if (uactive.getTiltedParticle() != null)
					p = uactive.getTiltedParticle();
				else
				{
					p = new TiltedParticle(x, y, uc.getActiveParticle());

					uc.getActiveParticle().setTiltedParticle(p);
					um.getTiltedMicrograph().addParticle(p);
				}
				active = p;
				frame.updateMicrographsModel();
			}
			frame.setChanged(true);
			repaint();
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
		if (SwingUtilities.isRightMouseButton(e))
			return;//dont move particle if it is dragging image
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (active != null && um.fits(x, y, frame.getParticleSize()))
		{
			setActiveMoved(true);
			moveActiveParticle(x, y);
			
		}
		frame.setChanged(true);
		repaint();
	}

	public void mouseReleased(MouseEvent e)
	{
		super.mouseReleased(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		manageActive(x, y);
	}


	@Override
	public void refreshActive(Particle p)
	{
		if(p != null)
			frame.getCanvas().refreshActive(((TiltedParticle) p).getUntiltedParticle());
	}

	@Override
	public ParticlePickerJFrame getFrame()
	{
		return frame;
	}

	@Override
	public Micrograph getMicrograph()
	{
		return um.getTiltedMicrograph();
	}

	@Override
	public TrainingParticle getActive()
	{
		return active;
	}

	@Override
	protected void doCustomPaint(Graphics2D g2)
	{
		g2.setColor(frame.getColor());
		int index = 0;
		List<TiltedParticle> particles = um.getTiltedMicrograph().getParticles();
		for (TiltedParticle p : particles)
		{
			drawShape(g2, p, index == (particles.size() - 1));
			index++;
		}
		
		if (uc.getActiveTiltedParticle() != null)
		{
			g2.setColor(Color.red);
			drawShape(g2, uc.getActiveParticle().getTiltedParticle(), true);
		}
		if(frame.drawAngles())
			drawLine(Math.toRadians(um.getTiltedAngle()), g2);
		
	}

	@Override
	protected Particle getLastParticle()
	{
		return uc.getLastParticle().getTiltedParticle();
	}

	
	protected void manageActive(int x, int y)
	{
		if(!activemoved)
			return;
		if (active.getUntiltedParticle().isAdded())
			um.initAligner();
		setActiveMoved(false);
	}
	
	@Override
	public void setMicrograph(Micrograph m)
	{
		um = (UntiltedMicrograph)m;
	}
	
	


}
