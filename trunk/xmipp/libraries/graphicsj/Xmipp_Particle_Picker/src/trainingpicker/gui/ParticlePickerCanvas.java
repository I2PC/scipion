package trainingpicker.gui;

import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.BasicStroke;
import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.List;

import javax.swing.SwingUtilities;

import trainingpicker.model.AutomaticParticle;
import trainingpicker.model.FamilyState;
import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.MicrographFamilyData;
import trainingpicker.model.MicrographFamilyState;
import trainingpicker.model.TrainingParticle;
import trainingpicker.model.TrainingPicker;
import xmipp.Particle;

public class ParticlePickerCanvas extends ImageCanvas implements MouseWheelListener
{

	private ParticlePickerJFrame frame;
	private TrainingMicrograph micrograph;
	private TrainingParticle dragged;
	private TrainingPicker ppicker;
	final static BasicStroke dashedst = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 10.0f }, 0.0f);
	final static BasicStroke continuousst = new BasicStroke();

	public ParticlePickerCanvas(ParticlePickerJFrame frame)
	{
		super(frame.getMicrograph().getImage());
		this.micrograph = frame.getMicrograph();
		this.frame = frame;
		addMouseWheelListener(this);
		this.ppicker = frame.getParticlePicker();

	}

	public void updateMicrograph()
	{
		this.micrograph = frame.getMicrograph();
		ImageWindow iw = (ImageWindow) getParent();
		iw.setImage(micrograph.getImage());
		iw.updateImage(micrograph.getImage());
		iw.setTitle(micrograph.getName());
		imp = micrograph.getImage();
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
			return;
		}
		if (frame.getFamilyData().isPickingAvailable())
		{
			TrainingParticle p = micrograph.getParticle(x, y);
			if (p == null)
				p = micrograph.getAutomaticParticle(x, y, frame.getThreshold());
			if (p != null)
			{
				if (SwingUtilities.isLeftMouseButton(e) && e.isControlDown())
				{
					micrograph.removeParticle(p, ppicker);
					frame.updateMicrographsModel();
				}
				else if (SwingUtilities.isLeftMouseButton(e))
				{
					dragged = p;
				}
			}
			else if (SwingUtilities.isLeftMouseButton(e) && TrainingParticle.boxContainedOnImage(x, y, frame.getFamily().getSize(), imp))
			{
				p = new TrainingParticle(x, y, frame.getFamily(), micrograph);
				micrograph.addManualParticle(p);
				dragged = p;
				frame.updateMicrographsModel();
			}
			frame.setChanged(true);
			repaint();
		}
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
		dragged = null;
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
			return;
		}
		if (frame.getFamilyData().isPickingAvailable())
		{
			if (dragged == null)
				return;

			if (!TrainingParticle.boxContainedOnImage(x, y, dragged.getFamily().getSize(), imp))
				return;
			if (dragged instanceof AutomaticParticle)
			{
				micrograph.removeParticle(dragged, ppicker);
				dragged = new TrainingParticle(dragged.getX(), dragged.getY(), dragged.getFamily(), micrograph);
				micrograph.addManualParticle(dragged);
			}
			else
			{
				dragged.setPosition(x, y);
				if (frame.getParticlesJDialog() != null)
					dragged.getParticleCanvas(frame).repaint();
			}
			frame.setChanged(true);
			repaint();
		}
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);

	}

	public void paint(Graphics g)
	{
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		if (frame.getFamily().getStep() == FamilyState.Manual)
			for (MicrographFamilyData mfdata : micrograph.getFamiliesData())
				drawFamily(g2, mfdata);
		else
			drawFamily(g2, micrograph.getFamilyData(frame.getFamily()));
		Rectangle r = getSrcRect();
        g2.drawRect(0, 0, (int)(r.getWidth() * magnification - 1), (int)(r.getHeight() * magnification - 1));					
	}

	private void drawFamily(Graphics2D g2, MicrographFamilyData mfdata)
	{
		int x0 = (int) getSrcRect().getX();
		int y0 = (int) getSrcRect().getY();

		int radius;
		List<TrainingParticle> particles;
		int index;
		if (!mfdata.isEmpty())
		{
			particles = mfdata.getManualParticles();
			g2.setColor(mfdata.getFamily().getColor());
			radius = (int) (mfdata.getFamily().getSize() / 2 * magnification);

			for (index = 0; index < particles.size(); index++)
				drawShape(g2, particles.get(index), x0, y0, radius, index == particles.size() - 1);
			List<AutomaticParticle> autoparticles = mfdata.getAutomaticParticles();
			for (int i = 0; i < autoparticles.size(); i++)
				if (!autoparticles.get(i).isDeleted() && autoparticles.get(i).getCost() >= frame.getThreshold())
					drawShape(g2, autoparticles.get(i), x0, y0, radius, false);
		}
	}

	private void drawShape(Graphics2D g2, TrainingParticle p, int x0, int y0, int radius, boolean all)
	{

		int x = (int) ((p.getX() - x0) * magnification);
		int y = (int) ((p.getY() - y0) * magnification);
		int distance = (int) (10 * magnification);
		if (p instanceof AutomaticParticle)
			g2.setStroke(dashedst);
		if (frame.isShapeSelected(Shape.Rectangle) || all)
			g2.drawRect(x - radius, y - radius, radius * 2, radius * 2);
		if (frame.isShapeSelected(Shape.Circle) || all)
			g2.drawOval(x - radius, y - radius, radius * 2, radius * 2);
		g2.setStroke(continuousst);
		if (frame.isShapeSelected(Shape.Center) || all)
		{
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}
	}

	public void moveTo(TrainingParticle p)
	{
		int width = (int) getSrcRect().getWidth();
		int height = (int) getSrcRect().getHeight();
		int x0 = p.getX() - width / 2;
		if(x0 < 0)
			x0 = 0;
		if(x0 + width > imp.getWidth())
			x0 = imp.getWidth() - width;
		int y0 = p.getY() - height / 2;
		if(y0 < 0)
			y0 = 0;
		if(y0 + height > imp.getHeight())
			y0 = imp.getHeight() - height;
		Rectangle r = new Rectangle(x0, y0, width, height);
		if (!getSrcRect().contains(r))
		{
			setSourceRect(r);
			repaint();
		}
	}

}
