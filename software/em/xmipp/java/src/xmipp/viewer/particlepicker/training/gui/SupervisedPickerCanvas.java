package xmipp.viewer.particlepicker.training.gui;

import java.awt.*;
import java.awt.event.KeyEvent;
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

    public static final int POLIGONAL_MODE_KEY = KeyEvent.VK_C;
    public static final int ERASER_MEDIUM_SIZE_KEY = KeyEvent.VK_M;
    public static final int ERASER_LARGE_SIZE_KEY = KeyEvent.VK_B;
    public static final int ERASER_X_LARGE_SIZE_KEY = KeyEvent.VK_X;
    private ManualParticle active;
    private Point linearPickingStartPoint;


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
				erase(x, y, e);
				return;
			}

            // Filament/linear picking
            if (frame.isLinearMode()) {
                linearPickingMousePressed(x,y, isKeyPressed(POLIGONAL_MODE_KEY));
                return;
            };

			ManualParticle p = getMicrograph().getParticle(x, y);
			if (p == null)
				p = getMicrograph().getAutomaticParticle(x, y, getFrame().getThreshold());

            if (p != null && SwingUtilities.isLeftMouseButton(e))
			{
				active = p;
				repaint();

			}
			else if (SwingUtilities.isLeftMouseButton(e)){
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

            if(getFrame().isCenterParticle())
                new CenterParticleTask(this, getParticlePicker(), p).execute();

            active = p;
            if(picker.getMode() == Mode.Manual)
                new ParticleToTemplatesTask(active).execute();

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

    protected void erase(int x, int y)
	{
        erase(x,y,true);
	}

    protected void erase(int x, int y, boolean refresh)
    {
        getMicrograph().removeParticles(x, y, getParticlePicker());
        active = getLastParticle();
        if (refresh) refresh();
    }

    // Erases having a larger point (square in this case)
    protected void erase(int x, int y, MouseEvent e)
    {

        int particleSize = getFrame().getSide();

        //
        int size = 1;

        if (isKeyPressed(ERASER_MEDIUM_SIZE_KEY)) {
            size = 2;
        } else if (isKeyPressed(ERASER_LARGE_SIZE_KEY)){
            size = 3;
        } else if (isKeyPressed(ERASER_X_LARGE_SIZE_KEY)){
            size = 4;
        }

        // Calculate the eraser "area" (square)
        int side = (int) Math.ceil((size*particleSize/10)/getMagnification());

        // Make it odd
        if (side % 2 == 0){
            side = side + 1;
        }

        // Half of if
        int half = (side -1)/2;

        int step = Math.round(particleSize/2);

        for (int x2 = x-half; x2<=x+half; x2+=step ){

            for (int y2 = y-half; y2<=y+half; y2+=step ){
                erase(x2,y2, false);
            }
        }

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
					erase(x, y, e);
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
				erase(x, y, e);
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

    public void mouseMoved(MouseEvent e) {
        super.mouseMoved(e);

        if (getFrame().isLinearMode()){
            paintLine(e);
        }
    }

    public void keyPressed(KeyEvent e){

        super.keyPressed(e);

        if (e.getKeyCode() == KeyEvent.VK_ESCAPE){

            if (haveFirstPoint()){
                resetLinearPicking();
            }
        }
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
