package xmipp.viewer.particlepicker.training.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.util.List;
import javax.swing.SwingUtilities;
import xmipp.jni.Particle;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.training.model.AutomaticParticle;
import xmipp.viewer.particlepicker.training.model.CenterParticleTask;
import xmipp.viewer.particlepicker.training.model.ManualParticle;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.ParticleToTemplatesTask;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

public class SupervisedPickerCanvas extends ParticlePickerCanvas
{

	
	private ManualParticle active;

	public SupervisedPickerCanvas(SupervisedPickerJFrame frame)
	{
		super(frame);
		
		active = getLastParticle();

	}

	protected ManualParticle getLastParticle()

	{
		if (getMicrograph().getParticles().isEmpty())
			return null;
		return getMicrograph().getLastAvailableParticle(getFrame().getThreshold());
	}

	public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);
		int x = offScreenX(e.getX());
		int y = offScreenY(e.getY());

		if (frame.isPickingAvailable(e))
		{
			if (frame.isEraserMode())
			{
				erase(x, y);
				return;
			}

			ManualParticle p = getMicrograph().getParticle(x, y);
			if (p == null)
				p = getMicrograph().getAutomaticParticle(x, y, getFrame().getThreshold());
			if (p != null && SwingUtilities.isLeftMouseButton(e))
			{
				active = p;
				repaint();

			}
			else if (SwingUtilities.isLeftMouseButton(e) && micrograph.fits(x, y, frame.getParticlePicker().getSize()))
			{
                            
				p = new ManualParticle(x, y, frame.getParticlePicker(), micrograph);
				getMicrograph().addManualParticle(p, getParticlePicker());
                                if(getFrame().isCenterParticle())
                                    new CenterParticleTask(this, getParticlePicker(), p).execute();
				active = p;
				if(picker.getMode() == Mode.Manual)
//					TasksManager.getInstance().addTask(new ParticleToTemplatesTask(active));
					new ParticleToTemplatesTask(active).execute();
				refresh();
			}
		}
	}

	protected void erase(int x, int y)
	{
		getMicrograph().removeParticles(x, y, getParticlePicker());
		active = getLastParticle();
		refresh();
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	@Override
	public void mouseDragged(MouseEvent e)
	{

		super.mouseDragged(e);
		if (SwingUtilities.isLeftMouseButton(e))
		{
			int x = offScreenX(e.getX());
			int y = offScreenY(e.getY());
			if (frame.isPickingAvailable(e))
			{
				if (frame.isEraserMode())
				{
					erase(x, y);
					return;
				}
				if (active == null)
					return;

				if (!micrograph.fits(x, y, active.getParticlePicker().getSize()))
					return;
				if (active instanceof AutomaticParticle)
				{
					getMicrograph().removeParticle(active, getParticlePicker());
					active = new ManualParticle(active.getX(), active.getY(), picker, micrograph);
					getMicrograph().addManualParticle(active, getParticlePicker());
					
				}
				else
				{
					setActiveMoved(true);
					moveActiveParticle(x, y);
				}

				manageActive(x, y);

			}
		}
	}

	@Override
	public void mouseReleased(MouseEvent e)
	{

		super.mouseReleased(e);
		int x = offScreenX(e.getX());
		int y = offScreenY(e.getY());
		if (frame.isPickingAvailable(e))
		{
			if (frame.isEraserMode())
			{
				erase(x, y);
				return;
			}
			//deleting when mouse released, takes less updates to templates and frame
			if (active != null && SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
			{
				getMicrograph().removeParticle(active, getParticlePicker());
				active = getLastParticle();
				refresh();

			}
			else
				manageActive(x, y);
			if (activemoved)
			{
				setActiveMoved(false);
				if(picker.getMode() == Mode.Manual)
//					TasksManager.getInstance().addTask(new ParticleToTemplatesTask(active));
					new ParticleToTemplatesTask(active).execute();
			}
		}
	}

	public void manageActive(int x, int y)
	{
		if (!activemoved)
			return;

		if (!micrograph.fits(x, y, active.getSize()))
			return;

		if (active instanceof AutomaticParticle && !((AutomaticParticle) active).isDeleted())
		{

			getMicrograph().removeParticle(active, getParticlePicker());
			active = new ManualParticle(active.getX(), active.getY(), active.getParticlePicker(), micrograph);
			getMicrograph().addManualParticle(active, getParticlePicker());
			repaint();
		}
		else
		{
			//changing particle position requires to remove it from template and add it again
			moveActiveParticle(x, y);
			repaint();

		}
		frame.setChanged(true);
	}

	protected void doCustomPaint(Graphics2D g2)
	{
		List<ManualParticle> particles;
		int index;
		Color color = picker.getColor();
		Color autoColor = color.darker();
		if (!getMicrograph().isEmpty())
		{
			particles = getMicrograph().getManualParticles();
			g2.setColor(color);

			for (index = 0; index < particles.size(); index++)
				drawShape(g2, particles.get(index), false, continuousst);

			g2.setColor(autoColor);
			List<AutomaticParticle> autoparticles = getMicrograph().getAutomaticParticles();
			for (int i = 0; i < autoparticles.size(); i++)
				if (!autoparticles.get(i).isDeleted() && autoparticles.get(i).getCost() >= getFrame().getThreshold())
					drawShape(g2, autoparticles.get(i), false, continuousst);

		}
		if (active != null)
		{
			boolean isauto = active instanceof AutomaticParticle;
			
			color = isauto? Color.red.darker(): Color.red;
			g2.setColor(color);
			drawShape(g2, active, true, activest);
		}
		Rectangle autopickout = getMicrograph().getRectangle();
		if (autopickout != null && getMicrograph().hasManualParticles())
		{
			g2.setColor(Color.yellow);
			g2.setStroke(continuousst);
			int x = getXOnImage((int) autopickout.getX());
			int y = getYOnImage((int) autopickout.getY());
			g2.drawRect(x, y, (int) (autopickout.getWidth() * magnification), (int) (autopickout.getHeight() * magnification));
		}

	}

	@Override
	public void refreshActive(Particle p)
	{
		if (p == null)
			active = null;
		else
			active = (ManualParticle) p;
		repaint();

	}

	@Override
	public SupervisedPickerJFrame getFrame()
	{
		return (SupervisedPickerJFrame)frame;
	}

	@Override
	public SupervisedPickerMicrograph getMicrograph()
	{
		return (SupervisedPickerMicrograph)micrograph;
	}

	@Override
	public ManualParticle getActive()
	{
		return active;
	}

	
	
	@Override
	public SupervisedParticlePicker getParticlePicker()
	{
		// TODO Auto-generated method stub
		return (SupervisedParticlePicker)picker;
	}

}
