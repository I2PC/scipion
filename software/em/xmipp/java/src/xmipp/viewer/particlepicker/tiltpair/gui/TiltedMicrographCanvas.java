package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.*;
import java.awt.event.MouseEvent;

import java.util.List;

import javax.swing.SwingUtilities;

import xmipp.jni.Particle;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.*;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;

public class TiltedMicrographCanvas extends ParticlePickerCanvas
{

    // Fields
    // Untilted references
    private UntiltedMicrograph um;
    private UntiltedMicrographCanvas uc;

    // Constructors
    public TiltedMicrographCanvas(TiltPairPickerJFrame frame)
    {
        super(frame, frame.getMicrograph().getTiltedMicrograph());
        this.um = frame.getMicrograph();
        this.uc = (UntiltedMicrographCanvas) frame.getCanvas();

    }

    // Fields getters and setters
    @Override
    public TiltPairPickerJFrame getFrame() {return (TiltPairPickerJFrame)frame;}

    @Override
    public TiltedParticle getActive() { return uc.getActiveTiltedParticle();}

    @Override
    public TiltedMicrograph getMicrograph() {return (TiltedMicrograph)micrograph;}

    @Override
    protected TiltedParticle getLastParticle()
    {
        return uc.getLastParticle().getTiltedParticle();
    }

    @Override
    public TiltPairPicker getParticlePicker()
    {
        return (TiltPairPicker)picker;
    }

    // Methods
	public void updateMicrograph()
	{
		um = getFrame().getMicrograph();
		TiltedMicrograph m = getFrame().getMicrograph().getTiltedMicrograph();
		updateMicrograph(m);
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

    private void addParticle(int x, int y)
    {
        try
        {

            Particle p = um.getAlignedUntiltedParticle(x, y);
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
    @Override
    protected void removeParticle(PickerParticle particleToRemove) {

        TiltedParticle toRemove = (TiltedParticle) particleToRemove;
        uc.removeParticle(toRemove.getUntiltedParticle());

    }

    @Override
    public void refreshActive(PickerParticle p)
    {
        if (p != null){
            active = p;
            frame.getCanvas().refreshActive(((TiltedParticle) p).getUntiltedParticle());
        }
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
            drawShape(g2, uc.getActive().getTiltedParticle(), true);
        }
        if (getFrame().drawAngles())
            drawLine(Math.toRadians(um.getTiltedAngle()), g2);

    }

    // Event listeners

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
            if (frame.isEraserMode())
            {
                eraseAndSync(e);
                return;
            }

            int x = super.offScreenX(e.getX());
            int y = super.offScreenY(e.getY());

            TiltedParticle p = um.getTiltedMicrograph().getParticle(x, y, getFrame().getParticleSize());
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
                    UntiltedParticle uactive = uc.getActive();
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

    private void eraseAndSync(MouseEvent e) {
        erase(e);
        getFrame().getCanvas().repaint();

    }

    public void mousePressed(int x, int y)
	{
		setupScroll(x, y);
	}

	/**
	 * Syncs tilted micrograph with untilted
	 */
	public void syncWithUntiltedCanvas(int x, int y)
	{

        scroll(x,y);

	}

	@Override
	public void mouseDragged(MouseEvent e)
	{
		super.mouseDragged(e);

		if (frame.isPickingAvailable(e)){

            if (frame.isEraserMode())
            {
                eraseAndSync(e);
                return;
            }

            int x = super.offScreenX(e.getX());
            int y = super.offScreenY(e.getY());

            if (uc.getActive().getTiltedParticle() != null && um.fits(x, y, getFrame().getParticleSize())){
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

    public void setupScrollFromUntilted(int sx, int sy) {

        // Calculate the offset coordinates for sx and sy
        int ox = offScreenX(sx);
        int oy = offScreenY(sy);

        setupScroll(ox, oy);

    }
}
