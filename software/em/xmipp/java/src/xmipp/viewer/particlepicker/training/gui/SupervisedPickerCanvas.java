package xmipp.viewer.particlepicker.training.gui;

import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.util.List;
import javax.swing.SwingUtilities;


import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.PickerParticle;
import xmipp.viewer.particlepicker.training.model.AutomaticParticle;
import xmipp.viewer.particlepicker.training.model.CenterParticleTask;
import xmipp.viewer.particlepicker.training.model.ManualParticle;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.ParticleToTemplatesTask;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

public class SupervisedPickerCanvas extends ParticlePickerCanvas
{

    public static final int POLIGONAL_MODE_KEY = KeyEvent.VK_C;
    private Point linearPickingStartPoint;


    public SupervisedPickerCanvas(SupervisedPickerJFrame frame)
	{
		super(frame);
		
		active = getLastManualParticle();

	}

	protected ManualParticle getLastManualParticle(){

        if (getMicrograph().getParticles().isEmpty())
			return null;
		return getMicrograph().getLastAvailableParticle(getFrame().getThreshold());
	}
    @Override
    public ManualParticle getActive(){
        return (ManualParticle) active;
    }

    protected ManualParticle getLastParticle(){
        return getLastManualParticle();
    }


    public void mousePressed(MouseEvent e)
	{
		super.mousePressed(e);

		int x = offScreenX(e.getX());
		int y = offScreenY(e.getY());

		if (frame.isPickingAvailable(e))
		{
			if (frame.isEraserMode(e)){
				erase(e);
				return;
			}

            // Filament/linear picking
            if (frame.isLinearMode()) {
                linearPickingMousePressed(x,y, isKeyPressed(POLIGONAL_MODE_KEY));
                return;
            }

			ManualParticle p = getMicrograph().getParticle(x, y);
			if (p == null)
				p = getMicrograph().getAutomaticParticle(x, y, getFrame().getThreshold());

            if (p != null && SwingUtilities.isLeftMouseButton(e)){
                active = p;
				repaint();

			}else if (SwingUtilities.isLeftMouseButton(e)){
                addParticle(x, y, true);
            }
        }
	}

    protected void addParticle( int x, int y, boolean refresh) {

        ManualParticle p;
        if (micrograph.fits(x, y, frame.getParticlePicker().getSize()))
        {

            p = new ManualParticle(x, y, frame.getParticlePicker(), micrograph);
            getMicrograph().addManualParticle(p, getParticlePicker());

            if(getFrame().centerParticle())
                new CenterParticleTask(this, getParticlePicker(), p).execute();

            active = p;
            if(picker.getMode() == Mode.Manual)
                new ParticleToTemplatesTask(getActive()).execute();

            if (refresh) refresh();
        }
    }

    private void linearPickingMousePressed(int x, int y, boolean continuous){

        Point point = new Point(x,y);

        if (haveFirstPoint()) {
            finishLinearPicking(point, continuous);
        } else {
            startLinearPicking(point);
        }
    }

    private boolean haveFirstPoint() {
        return linearPickingStartPoint != null;
    }

    private void startLinearPicking(Point point) {

        this.linearPickingStartPoint = point;
    }

    private void finishLinearPicking(Point point, boolean continuous) {

        addParticlesInLine(point);

        // Start from here.
        if (continuous) {
            this.linearPickingStartPoint = point;
        } else {
            resetLinearPicking();
        }

    }

    public void resetLinearPicking(){
        this.linearPickingStartPoint = null;
    }

    private void addParticlesInLine(Point endPoint) {

        // To be moved to the interface!!
        int gap = this.getFrame().getStep();

        // Get the distance diff (x,y) between ent and start
        double lengthX = endPoint.getX() - this.linearPickingStartPoint.getX();
        double lengthY = endPoint.getY() - this.linearPickingStartPoint.getY();

        // Get the angle (tan-1: inverse)
        double angle = Math.atan2(lengthY, lengthX);

        // Calculate the number length of the line
        double lineLength = Math.sqrt(Math.pow(lengthY,2) + Math.pow(lengthX,2));

        // Calculate number of segments for the iteration
        int segmentsNum = (int) Math.round(lineLength/gap);

        double actualGap = lineLength/segmentsNum;

        // Now create as one particle each segment of size "gap"
        for (int step = 0; step <= segmentsNum; step++){

            // Calculate the X
            int x = (int)Math.round(this.linearPickingStartPoint.getX() + (Math.cos(angle)*step*actualGap));
            int y = (int)Math.round(this.linearPickingStartPoint.getY() + (Math.sin(angle)*step*actualGap));

            addParticle(x,y,step==segmentsNum);

        }

    }

	/**
	 * Updates particle position and repaints if on pick.
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
				if (frame.isEraserMode(e))
				{
					erase(e);
					return;
				}
				if (active == null)
					return;

                ManualParticle activeParticle = getActive();

				if (!micrograph.fits(x, y, activeParticle.getParticlePicker().getSize()))
					return;
				if (active instanceof AutomaticParticle)
				{
					getMicrograph().removeParticle(activeParticle, getParticlePicker());
					active = new ManualParticle(activeParticle.getX(), activeParticle.getY(), picker, micrograph);
					getMicrograph().addManualParticle((ManualParticle) active, getParticlePicker());
					refresh();
					
				}else{
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
			if (frame.isEraserMode(e))
			{
				erase(e);
				return;
			}
			//deleting when mouse released, takes less updates to templates and frame
			if (active != null && SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
			{
				getMicrograph().removeParticle(getActive(), getParticlePicker());
				active = getLastManualParticle();
				refresh();

			}
			else
				manageActive(x, y);
			if (activemoved)
			{
				setActiveMoved(false);
				if(picker.getMode() == Mode.Manual)
//					TasksManager.getInstance().addTask(new ParticleToTemplatesTask(active));
					new ParticleToTemplatesTask(getActive()).execute();
			}
		}
	}

	public void manageActive(int x, int y)
	{

        ManualParticle activeParticle = getActive();

		if (!activemoved)
			return;

		if (!micrograph.fits(x, y, activeParticle.getSize()))
			return;

		if (active instanceof AutomaticParticle && !((AutomaticParticle) active).isDeleted())
		{

			getMicrograph().removeParticle(activeParticle, getParticlePicker());
			active = new ManualParticle(activeParticle.getX(), activeParticle.getY(), activeParticle.getParticlePicker(), micrograph);
			getMicrograph().addManualParticle(activeParticle, getParticlePicker());
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
		SupervisedParticlePicker picker = getParticlePicker();
		Color color = picker.getColor();
        Color autoColor = picker.getAutomaticColor();
        Color delColor = picker.getDeletedColor();
        if (!getMicrograph().isEmpty()) {
            particles = getMicrograph().getManualParticles();
            g2.setColor(color);

            for (index = 0; index < particles.size(); index++)
                drawShape(g2, particles.get(index), false, thinContinuousSt);

            List<AutomaticParticle> autoparticles = getMicrograph().getAutomaticParticles();
            for (AutomaticParticle autoparticle : autoparticles)
                if (!autoparticle.isDeleted() && autoparticle.getCost() >= getFrame().getThreshold()){
                    g2.setColor(autoColor);
                    drawShape(g2, autoparticle, false, thinContinuousSt);
                }else if (autoparticle.isUnavailable() && picker.isShowDeleted()) {
                    g2.setColor(delColor);
                    drawShape(g2, autoparticle, false, thinContinuousSt);
                }
		}
		if (active != null)
		{
			boolean isauto = active instanceof AutomaticParticle;
			
			color = isauto? Color.red.darker(): Color.red;
			g2.setColor(color);
			drawShape(g2, getActive(), true, activest);
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
	public void refreshActive(PickerParticle p)
	{
		if (p == null)
			// active = null;
            active = getLastParticle();
		else
			active = p;
		repaint();

	}

	@Override
	public SupervisedPickerJFrame getFrame()
	{
		return (SupervisedPickerJFrame)frame;
	}

    @Override
    protected void removeParticle(PickerParticle particleToRemove) {
        getMicrograph().removeParticle(particleToRemove,getParticlePicker());
    }

    @Override
	public SupervisedPickerMicrograph getMicrograph()
	{
		return (SupervisedPickerMicrograph)micrograph;
	}
	
	
	@Override
	public SupervisedParticlePicker getParticlePicker()
	{
		// TODO Auto-generated method stub
		return (SupervisedParticlePicker)picker;
	}


    public void keyPressed(KeyEvent e){

        super.keyPressed(e);

        if (e.getKeyCode() == KeyEvent.VK_ESCAPE){

            if (haveFirstPoint()){
                resetLinearPicking();
            }
        }
    }

    public void mouseMoved(MouseEvent e) {
        super.mouseMoved(e);

        if (getFrame().isLinearMode())
            paintLine(e);
    }

    private void paintLine(MouseEvent e){

        if (haveFirstPoint()){

            Graphics g = this.getGraphics();

            // Repaint
            this.paint(g);
            // Draw a line from the first point to the current mouse cursor point.

            g.setColor(Color.YELLOW);
            g.drawLine(getXOnImage(this.linearPickingStartPoint.x), getYOnImage(this.linearPickingStartPoint.y), e.getX(),e.getY());
        }

    }

}
