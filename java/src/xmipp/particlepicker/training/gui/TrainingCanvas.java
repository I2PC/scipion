package xmipp.particlepicker.training.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Stroke;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.util.List;

import javax.swing.SwingUtilities;

import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePickerCanvas;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.training.model.AutomaticParticle;
import xmipp.particlepicker.training.model.FamilyState;
import xmipp.particlepicker.training.model.MicrographFamilyData;
import xmipp.particlepicker.training.model.TrainingMicrograph;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.particlepicker.training.model.TrainingPicker;

public class TrainingCanvas extends ParticlePickerCanvas 
{

	private TrainingPickerJFrame frame;
	private TrainingMicrograph micrograph;
	private TrainingParticle active;
	private TrainingPicker ppicker;
	final static BasicStroke dashedst = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 10.0f }, 0.0f);
	final static BasicStroke continuousst = new BasicStroke();
	final static BasicStroke activedst = new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 10.0f }, 0.0f);
	final static BasicStroke activecst = new BasicStroke(2.0f);
	

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
			setActive(getLastParticle());
		else
			active = null;
		
	}
	
	TrainingParticle getLastParticle()
	{
		return frame.getFamilyData().getParticles().get(frame.getFamilyData().getParticles().size() - 1);
	}


	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		if (frame.isPickingAvailable(e))
		{

			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());
			TrainingParticle p = micrograph.getParticle(x, y);
			if (p == null)
				p = micrograph.getAutomaticParticle(x, y, frame.getThreshold());
			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
				{
					micrograph.removeParticle(p, ppicker);
					frame.updateMicrographsModel();
					active = frame.getFamilyData().getLastAvailableParticle(frame.getThreshold());
				}
				else if (SwingUtilities.isLeftMouseButton(e))
					active = p;
			}
			else if (SwingUtilities.isLeftMouseButton(e) && micrograph.fits(x, y, frame.getFamily().getSize()))
			{
				p = new TrainingParticle(x, y, frame.getFamily(), micrograph);
				micrograph.addManualParticle(p);
				active = p;
				frame.updateMicrographsModel();
			}
			frame.setChanged(true);
			repaint();
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
				moveActiveParticle(x, y);
				if(frame.templatesdialog != null)
					frame.loadTemplates();
			}
			frame.setChanged(true);
			repaint();
		}
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
	public void setActive(TrainingParticle p)
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
