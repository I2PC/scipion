
package xmipp.viewer.particlepicker.training.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Stroke;
import java.awt.event.MouseEvent;
import java.util.List;
import javax.swing.SwingUtilities;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.ParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.model.AutomaticParticle;
import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.MicrographFamilyData;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public class TrainingCanvas extends ParticlePickerCanvas
{

	private TrainingPickerJFrame frame;
	private TrainingMicrograph micrograph;
	private TrainingParticle active;
	private TrainingPicker ppicker;

	

	public TrainingCanvas(TrainingPickerJFrame frame)
	{
		super(frame.getMicrograph().getImagePlus(frame.getParticlePicker().getFilters()));
		
		this.micrograph = frame.getMicrograph();
		this.frame = frame;
		this.ppicker = frame.getParticlePicker();
		micrograph.runImageJFilters(ppicker.getFilters());
		if(!frame.getFamilyData().getParticles().isEmpty())
			active = getLastParticle();
		else
			active = null;

	}

	public void updateMicrograph()
	{
		this.micrograph = frame.getMicrograph();
		updateMicrographData();
		if(!frame.getFamilyData().getParticles().isEmpty())
			refreshActive(getLastParticle());
		else
			active = null;
		
	}
	
	TrainingParticle getLastParticle()
	{
		return frame.getFamilyData().getLastAvailableParticle(frame.getThreshold());
	}


	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (frame.isPickingAvailable(e))
		{

			TrainingParticle p = micrograph.getParticle(x, y);
			if (p == null)
				p = micrograph.getAutomaticParticle(x, y, frame.getThreshold());
			if (p != null && SwingUtilities.isLeftMouseButton(e))
			{
				active = p;
				repaint();

			}
			else if (SwingUtilities.isLeftMouseButton(e) && micrograph.fits(x, y, frame.getFamily().getSize()))
			{
				p = new TrainingParticle(x, y, frame.getFamily(), micrograph);
				micrograph.addManualParticle(p);
				ppicker.addParticleToTemplates(p, frame.isCenterPick());
				active = p;
				refresh();
			}
		}
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

		if (frame.isPickingAvailable(e))
		{
			if (active == null)
				return;

			if (!micrograph.fits(x, y, active.getFamily().getSize()))
				return;
			if (active instanceof AutomaticParticle)
			{
				micrograph.removeParticle(active, ppicker);
				active = new TrainingParticle(active.getX(), active.getY(), active.getFamily(), micrograph);
				micrograph.addManualParticle(active);
			}
			else
			{
				setActiveMoved(true);
				moveActiveParticle(x, y);
			}
			frame.setChanged(true);
			repaint();
		}
	}
	
	@Override
	public void mouseReleased(MouseEvent e)
	{

		super.mouseReleased(e);
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());
		if (frame.isPickingAvailable(e))
		{
			
			//deleting when mouse released, takes less updates to templates and frame
			if (active != null && SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
			{
				micrograph.removeParticle(active, ppicker);
				active = getLastParticle();
				refresh();
				if (!(active instanceof AutomaticParticle))
					frame.updateTemplates();
			}
			else
				manageActive(x, y);
			if (activemoved)
			{
				frame.updateTemplates();
				setActiveMoved(false);
			}

		}
	}
	
	public void manageActive(int x, int y)
	{
		if (!activemoved)
			return;

		if (!micrograph.fits(x, y, active.getFamily().getSize()))
			return;
		
		if (active instanceof AutomaticParticle && !((AutomaticParticle)active).isDeleted())
		{

			micrograph.removeParticle(active, ppicker);
			active = new TrainingParticle(active.getX(), active.getY(), active.getFamily(), micrograph);
			micrograph.addManualParticle(active);
			repaint();
		}
		else
		{
			moveActiveParticle(x, y);
			repaint();

		}
		frame.setChanged(true);
		setActiveMoved(false);
	}

	
	

	public void paint(Graphics g)
	{
		Graphics offgc;
		Image offscreen = null;
		Dimension d = getSize();

		// create the offscreen buffer and associated Graphics
		offscreen = createImage(d.width, d.height);
		offgc = offscreen.getGraphics();
		super.paint(offgc);
		Graphics2D g2 = (Graphics2D) offgc;
		if (frame.getFamily().getStep() == FamilyState.Manual)
			for (MicrographFamilyData mfdata : micrograph.getFamiliesData())
				drawFamily(g2, mfdata);
		else
			drawFamily(g2, micrograph.getFamilyData(frame.getFamily()));
		if (active != null)
		{
			g2.setColor(Color.red);
			BasicStroke activest = (active instanceof AutomaticParticle)? activedst: activecst;
			g2.setStroke(activest);
			drawShape(g2, active, true);
		}
		g.drawImage(offscreen, 0, 0, this);
	}

	private void drawFamily(Graphics2D g2, MicrographFamilyData mfdata)
	{
		List<TrainingParticle> particles;
		int index;
		if (!mfdata.isEmpty())
		{
			particles = mfdata.getManualParticles();
			g2.setColor(mfdata.getFamily().getColor());

			for (index = 0; index < particles.size(); index++)
				drawShape(g2, particles.get(index), index == particles.size() - 1);
			Stroke previous = g2.getStroke();
			g2.setStroke(dashedst);
			List<AutomaticParticle> autoparticles = mfdata.getAutomaticParticles();
			for (int i = 0; i < autoparticles.size(); i++)
				if (!autoparticles.get(i).isDeleted() && autoparticles.get(i).getCost() >= frame.getThreshold())
					drawShape(g2, autoparticles.get(i), false);
			g2.setStroke(previous);
		}
	}

	@Override
	public void refreshActive(TrainingParticle p)
	{
		active = p;
		repaint();

	}

	@Override
	public ParticlePickerJFrame getFrame()
	{
		return frame;
	}

	@Override
	public Micrograph getMicrograph()
	{
		return micrograph;
	}

	@Override
	public TrainingParticle getActive()
	{
		return active;
	}

	

}
