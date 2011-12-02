package particlepicker.tiltpair.gui;

import ij.gui.ImageWindow;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseWheelListener;
import java.util.List;
import javax.swing.SwingUtilities;

import particlepicker.Micrograph;
import particlepicker.ParticlePickerCanvas;
import particlepicker.ParticlePickerJFrame;
import particlepicker.WindowUtils;
import particlepicker.tiltpair.model.TiltedMicrograph;
import particlepicker.tiltpair.model.TiltedParticle;
import particlepicker.tiltpair.model.UntiltedMicrograph;
import particlepicker.tiltpair.model.UntiltedParticle;
import particlepicker.training.model.TrainingMicrograph;
import particlepicker.training.model.TrainingParticle;
import xmipp.Particle;

public class TiltedMicrographCanvas extends ParticlePickerCanvas
{

	private TiltPairPickerJFrame frame;
	private UntiltedMicrograph um;
	private UntiltedMicrographCanvas uc;
	private TiltedParticle active;
	private ImageWindow iw;
	private boolean reload;
	private boolean drawalpha;

	public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
	{
		super(frame.getMicrograph().getTiltedMicrograph().getImagePlus());
		this.um = frame.getMicrograph();
		this.frame = frame;
		this.uc = (UntiltedMicrographCanvas)frame.getCanvas();
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0.7, 0, iw);
		addMouseWheelListener(this);

	}

	public void updateMicrograph()
	{
		this.um = frame.getMicrograph();
		updateMicrographData();
		active = null;
	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	public void mousePressed(int x, int y)
	{
		setupScroll(x, y);
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
	public void mouseDragged(int x, int y)
	{
		scroll(x, y);
	}

	public void mouseWheelMoved(int x, int y, int rotation)
	{

		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		if (getMagnification() <= 1.0)
			imp.repaintWindow();

	}

	public void paint(Graphics g)
	{
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
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
			drawLine(frame.getTiltedAngle(), g2);
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
				if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
				{
					um.removeParticle(p.getUntiltedParticle());
					frame.updateMicrographsModel();
					if (p.getUntiltedParticle().isAdded())
						reload = true;
					frame.getCanvas().repaint();
				}
				else if (SwingUtilities.isLeftMouseButton(e))
					active = p;
			}
			else if (uc.hasActiveParticle() && SwingUtilities.isLeftMouseButton(e) && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
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

	@Override
	public void mouseDragged(MouseEvent e)
	{

		super.mouseDragged(e);
		if (SwingUtilities.isRightMouseButton(e))
			return;//dont move particle if it is dragging image
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (active != null && Particle.boxContainedOnImage(x, y, frame.getParticleSize(), imp))
		{
			active.setPosition(x, y);
			if (active.getUntiltedParticle().isAdded())
				reload = true;
		}
		frame.setChanged(true);
		repaint();
	}



	@Override
	public void setActive(TrainingParticle p)
	{
		frame.getCanvas().setActive(((TiltedParticle) p).getUntiltedParticle());
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
	


}
