package xmipp.viewer.particlepicker;

import ij.gui.ImageWindow;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.Toolkit;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.util.ArrayList;
import java.util.List;

import xmipp.ij.commons.XmippImageCanvas;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.Particle;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippResource;
import xmipp.viewer.particlepicker.training.model.ManualParticle;

import javax.swing.*;

public abstract class ParticlePickerCanvas<P extends PickerParticle> extends XmippImageCanvas
{
	public final static BasicStroke dashedst = new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 10.0f }, 0.0f);
	public final static BasicStroke continuousst = new BasicStroke(2.0f);
	public final static BasicStroke activest = new BasicStroke(3.0f);
	//public final static BasicStroke activedst = new BasicStroke(3.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 10.0f }, 0.0f);
        
	protected ParticlePickerJFrame frame;
	protected Micrograph micrograph;
	protected ParticlePicker picker;
	protected XmippImageWindow iw;
	protected Cursor eraserCursor = Toolkit.getDefaultToolkit().createCustomCursor(XmippResource.getIcon("eraserC.png").getImage(), new Point(3, 31), "Eraser");
    protected Cursor linearCursor = Toolkit.getDefaultToolkit().createCustomCursor(XmippResource.getIcon("linearPickingC.png").getImage(), new Point(0, 0), "Linear");
    private byte pickingCursorIndex = 1;
    protected Cursor pickingCursor = getPickingCursor(pickingCursorIndex);
    protected P active;

    Dimension size;
	public void setMicrograph(Micrograph m)
	{
		micrograph = m;

	}

	protected boolean activemoved;

	private static boolean tongleSetSelected = false;
	
	public ParticlePickerCanvas(ParticlePickerJFrame frame)
	{
		this(frame, frame.getMicrograph());
	}

	public ParticlePickerCanvas(ParticlePickerJFrame frame, Micrograph micrograph)
	{
		super(micrograph.getImagePlus(frame.getParticlePicker().getFilters()));
		this.frame = frame;
		this.picker = frame.getParticlePicker();
		this.micrograph = micrograph;
		micrograph.runImageJFilters(picker.getFilters());

        addKeyListener(this);
		addMouseMotionListener(this);

		this.iw = new XmippImageWindow(getImage(), this, null);
		getFrame().displayZoom(getMagnification());
		iw.setTitle(getMicrograph().getName());

	}
	
	public void keyPressed(KeyEvent e)
	{
		super.keyPressed(e);

        P active = getActive();
		int step = 1;
		int code = e.getKeyCode();
//		int x = screenX(imp.getWidth()/2);
//		int y = screenY(imp.getHeight()/2);
		if(active != null)
		{
			if (code == KeyEvent.VK_UP)
			{
				setActiveMoved(true);
				manageActive(active.getX(), active.getY() - step);
			}
			else if (code == KeyEvent.VK_DOWN)
			{
				setActiveMoved(true);
				manageActive(active.getX(), active.getY() + step);
			}
			else if (code == KeyEvent.VK_LEFT)
			{
				setActiveMoved(true);
				manageActive(active.getX() - step, active.getY());
			}
			else if (code == KeyEvent.VK_RIGHT)
			{
				setActiveMoved(true);
				manageActive(active.getX() + step, active.getY());
			}
		}
		if (code == KeyEvent.VK_SPACE)
		{
			getFrame().circlechb.setSelected(tongleSetSelected);
			getFrame().rectanglechb.setSelected(tongleSetSelected);
			tongleSetSelected = !tongleSetSelected;
		}else if (code == ParticlePickerJFrame.TOGGLE_ERASE_MODE_KEY){
            // Toggle erase mode.
            getFrame().activateEraseMode();
        }else if (code == ParticlePickerJFrame.TOGGLE_LINEAR_MODE_KEY){
            // Toggle linear mode
            getFrame().activateLinearMode();
        }else if (code == ParticlePickerJFrame.TOGGLE_NORMAL_MODE_KEY){
            // Toggle linear mode
            getFrame().activateNormalMode();

        }else if (code == KeyEvent.VK_PLUS){
            changePickingCursor((byte) (pickingCursorIndex+1));

        }else if (code == KeyEvent.VK_MINUS){
            changePickingCursor((byte) (pickingCursorIndex-1));

        }else
			return;// do not repaint if not needed
		repaint();
	}

    private void changePickingCursor(byte index) {
        Cursor newCursor = getPickingCursor(index);

        if (newCursor != null){
            pickingCursorIndex = index;
            pickingCursor = newCursor;
            setCustomCursor();
        }
    }

    private Cursor getPickingCursor(int index) {

        String pickingCursorImageName = "pickingC" + index + ".png";

        ImageIcon icon = XmippResource.getIcon(pickingCursorImageName);

        if (icon.getIconHeight() == -1) {
            // Fall back on crosshair cursor
            return crosshairCursor;
        } else {
            return Toolkit.getDefaultToolkit().createCustomCursor(icon.getImage(), new Point(2, 0), "Picking");
        }

    }

    protected void captureWindowSize(){
        size = this.imp.getWindow().getSize();
        this.imp.getWindow().setIgnoreRepaint(true);

    }
    protected void restoreWindowSize(){
        if (size != null) {
            this.imp.getWindow().setSize(size);
            this.imp.getWindow().setIgnoreRepaint(false);
        }
    }

    public void zoomIn(int sx, int sy)
	{

        captureWindowSize();

		super.zoomIn(sx, sy);
		getFrame().displayZoom(getMagnification());

        restoreWindowSize();
	}
	
	public void zoomOut(int sx, int sy)
	{
        captureWindowSize();

		super.zoomOut(sx, sy);
		getFrame().displayZoom(getMagnification());

        restoreWindowSize();
	}

	protected abstract P getLastParticle();

	protected abstract void manageActive(int x, int y);

    public boolean hasActiveParticle()
    {
        return active != null;
    }

	protected void setActiveMoved(boolean b)
	{
		activemoved = b;
	}

	protected void refresh()
	{
		getFrame().updateMicrographsModel();
		getFrame().setChanged(true);
		repaint();
	}

	public ImageWindow getIw()
	{
		return iw;
	}

	public void moveTo(PickerParticle p)
	{
		int width = (int) getSrcRect().getWidth();
		int height = (int) getSrcRect().getHeight();
		int x0 = p.getX() - width / 2;
		if (x0 < 0)
			x0 = 0;
		if (x0 + width > imp.getWidth())
			x0 = imp.getWidth() - width;
		int y0 = p.getY() - height / 2;
		if (y0 < 0)
			y0 = 0;
		if (y0 + height > imp.getHeight())
			y0 = imp.getHeight() - height;
		Rectangle r = new Rectangle(x0, y0, width, height);
		if (!getSrcRect().contains(r))
		{
			setSourceRect(r);
			repaint();
		}
	}

	public void mouseEntered(MouseEvent e)
	{
		super.mouseEntered(e);
        this.requestFocus();
		setCustomCursor(e);
	}

	public void mouseExited(MouseEvent e)
	{
		super.mouseExited(e);
		setCursor(new Cursor(Cursor.DEFAULT_CURSOR));

	}

    public void setCustomCursor(){
        setCustomCursor(null);
    }
	protected void setCustomCursor(MouseEvent e)
	{
        boolean proceed = true;

		if (e != null) {
            proceed =  getFrame().isPickingAvailable(e);
        }

        if (proceed) {
            if (getFrame().isEraserMode())
                setCursor(eraserCursor);

            else if (getFrame().isLinearMode())
                setCursor(linearCursor);

            else
                setCursor(pickingCursor);
        }
	}

	public abstract void refreshActive(PickerParticle p);
	

	public abstract P getActive();

	
	public abstract ParticlePickerJFrame getFrame();


	protected void drawShape(Graphics2D g2, ManualParticle p, boolean all, Stroke stroke)
	{
		drawShape(g2, p.getX(), p.getY(), picker.getSize(), all, stroke);
	}
	
	protected void drawShape(Graphics2D g2, ManualParticle p, boolean all)
	{
		drawShape(g2, p, all, continuousst);
	}
	
	protected void drawShape(Graphics2D g2, int x, int y, int size, boolean all)
	{
		drawShape(g2, x, y, size, all, continuousst);
				
	}
	protected void drawShape(Graphics2D g2, int x, int y, int size, boolean all, Stroke stroke, Color color)
	{
		g2.setColor(color);
		drawShape(g2, x, y, size, all, stroke);
	}
	protected void drawShape(Graphics2D g2, int x, int y, int size, boolean all, Stroke stroke)
	{

		g2.setStroke(stroke);
		int length = (int) (size * magnification);
		int radius = (int) (size / 2. * magnification);
		x = getXOnImage(x);
		y = getYOnImage(y);
		int distance = Math.min(10, (int) (radius/4. * magnification));

		if (getFrame().isShapeSelected(Shape.Rectangle) || all)
			g2.drawRect(x - radius, y - radius, length, length);
		if (getFrame().isShapeSelected(Shape.Circle) || all)
			g2.drawOval(x - radius, y - radius, length, length);
		if (getFrame().isShapeSelected(Shape.Center) || all)
		{
			g2.setStroke(activest);
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}
	}


	protected void drawLine(double alpha, Graphics2D g2)
	{
		int width = imp.getWidth();
		int height = imp.getHeight();
		double m;
		double x1, y1, x2, y2;
		if (alpha != Math.PI / 2)
		{
			m = Math.tan(alpha - Math.PI / 2);

			double y = height / 2.f;
			double x = y / m;

			if (Math.abs(x) > width / 2.f)// cuts in image sides
			{
				x1 = width;// on image
				y1 = getYCutOnImage(m, width / 2.f);
				x2 = 0;
				y2 = getYCutOnImage(m, -width / 2.f);
			}
			else
			// cuts in image top and bottom
			{
				y1 = 0;
				x1 = getXCutOnImage(m, height / 2.f);
				y2 = height;
				x2 = getXCutOnImage(m, -height / 2.f);
			}
		}
		else
		{
			x1 = 0;
			y1 = y2 = height / 2.f;
			x2 = width;

		}
		Color ccolor = g2.getColor();
		g2.setColor(Color.yellow);
		g2.drawLine((int) (x1 * magnification), (int) (y1 * magnification), (int) (x2 * magnification), (int) (y2 * magnification));
		g2.setColor(ccolor);
	}
	
	
	public int getXOnImage(int x)
	{
		int x0 = (int) getSrcRect().getX();
		return (int) ((x - x0) * magnification);
	}
	
	public int getYOnImage(int y)
	{
		int y0 = (int) getSrcRect().getY();
		return (int) ((y - y0) * magnification);
	}

	private double getYCutOnImage(double m, double x)
	{
		int height = imp.getHeight();
		return height / 2.f - m * x;
	}


	private double getXCutOnImage(double m, double y)
	{
		int width = imp.getWidth();
		return y / m + width / 2.f;
	}
	
	public void updateMicrograph()
	{
		Micrograph m = getFrame().getMicrograph();
		updateMicrograph(m);
	}

	public void updateMicrograph(Micrograph m)
	{
		getMicrograph().releaseImage();
		setMicrograph(m);
		imp = m.getImagePlus(picker.getFilters());
		m.runImageJFilters(picker.getFilters());
		refreshActive(null);
		Rectangle previousRect = (Rectangle)srcRect.clone();
		int previousDstWidth = dstWidth;
		int previousDstHeight = dstHeight;
		if(iw.isClosing())
			iw = new XmippImageWindow(getImage(), this, null);
		iw.setImage(imp);
		double zoom = getParticlePicker().getZoom();
		int imageWidth = imp.getWidth();
		int imageHeight = imp.getHeight();
		if(zoom != -1)
		{
			
			setDrawingSize(previousDstWidth, previousDstHeight);
			if (previousRect.width < imageWidth && previousRect.height < imageHeight)
				setSourceRect(previousRect);
			setMagnification(zoom);
			iw.pack();
		}
		else
			iw.maximize();
		getFrame().displayZoom(getMagnification());
		iw.setTitle(getMicrograph().getName());
		
	}

    protected void erase(MouseEvent e)
    {
        // Get the eraser size and coordinates and apply zoom
        int size = getFrame().getEraserSize()/2;
        double distance = (size);

        // Add the radio of the particles: we want to eras particle if eraser touches the circle of the particle
        distance = distance + (getParticlePicker().getSize()/2);


        int cursorX = offScreenX(e.getX());
        int cursorY = offScreenY (e.getY());

        List<PickerParticle> particlesToDelete = new ArrayList<PickerParticle>();
        // Go through the list of particles
        for (PickerParticle particle : getMicrograph().getParticleList()) {

            // if particle is within the distance (inside eraser)
            int length = length(particle.getX(),particle.getY(), cursorX, cursorY);

            if (length <= distance){
                particlesToDelete.add(particle);
            }
        }

        // For each particle to delete
        for (PickerParticle particle : particlesToDelete) {
            removeParticle((P) particle);
        }

        active = getLastParticle();
        refresh();
        paintEraser(e);
    }

    private int length (int x, int y, int x1, int y1) {

        return (int) Math.ceil(Math.sqrt(Math.pow((x - x1), 2) + Math.pow((y - y1), 2)));

    }

    protected abstract void removeParticle(P particleToRemove);

    public abstract Micrograph getMicrograph();

	protected void moveActiveParticle(int x, int y)
	{
		if (active == null)
			return;
		try
		{
			active.setPosition(x, y);
		}
		catch (Exception e)
		{
			XmippDialog.showError(getFrame(), e.getMessage());
		}
		if (getFrame().getParticlesJDialog() != null)
			active.getParticleCanvas(getFrame()).repaint();
	}

	public void paint(Graphics g) {
		Graphics offgc;//Off screen graphic
		Image offscreen = null;//Off screen image
		Dimension d = getSize();

		// create the offscreen buffer and associated Graphics
		offscreen = createImage(d.width, d.height);
		offgc = offscreen.getGraphics();
		
		super.paint(offgc);//super paint in off screen
		//my paint in offscreen
		Graphics2D g2 = (Graphics2D) offgc;
		doCustomPaint(g2);
		//drawing offscreen image
		g.drawImage(offscreen, 0, 0, this);
	}

    public void mouseMoved(MouseEvent e) {
        super.mouseMoved(e);

        if (getFrame().isEraserMode()){
            paintEraser(e);
        }
    }

	protected abstract void doCustomPaint(Graphics2D g2);
         
    public abstract ParticlePicker getParticlePicker();

    protected void paintEraser(MouseEvent e){

        paintEraser(e.getX(),e.getY());

    }

    private void paintEraser(int x, int y){

        Graphics2D g = (Graphics2D) this.getGraphics();

        // Repaint
        this.paint(g);
        // Draw a circle from the center
        g.setColor(Color.CYAN);
        int r = getFrame().getEraserSize();

        // Apply magnification
        r = (int) Math.round((getMagnification()*r));

        int ovalX = x-(r/2);
        int ovalY = y-(r/2);

        g.setStroke(continuousst);//2 pixel width
        g.drawOval(ovalX,ovalY,r,r);

    }

    @Override
    public void mouseWheelMoved(MouseWheelEvent e) {
        // Resize the eraser brush
        if (!e.isShiftDown()){
            // if erase mode
            if (getFrame().isEraserMode()){
                getFrame().setEraserSize(getFrame().getEraserSize()+(e.getWheelRotation()*50));
                paintEraser(e);
            }

        } else {
            super.mouseWheelMoved(e);
        }
    }

}
