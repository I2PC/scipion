package xmipp.viewer.particlepicker.tiltpair.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;

import javax.swing.SwingUtilities;

import xmipp.jni.Particle;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippMessageDialog;
import xmipp.viewer.particlepicker.ParticlePickerCanvas;
import xmipp.viewer.particlepicker.PickerParticle;
import xmipp.viewer.particlepicker.tiltpair.model.TiltPairPicker;
import xmipp.viewer.particlepicker.tiltpair.model.TiltedParticle;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.viewer.particlepicker.training.model.ManualParticle;

public class UntiltedMicrographCanvas extends ParticlePickerCanvas
{
    // Constructors
    public UntiltedMicrographCanvas(TiltPairPickerJFrame frame)
    {
        super(frame);

    }

    // Fields getters and setters
	@Override
	public TiltPairPickerJFrame getFrame()
	{
		return (TiltPairPickerJFrame)frame;
	}

    @Override
	public UntiltedParticle getActive()
	{
		return (UntiltedParticle) active;
	}

    @Override
    public UntiltedMicrograph getMicrograph()
    {
        return (UntiltedMicrograph)micrograph;
    }

    protected UntiltedParticle getLastParticle()
    {
        if (getMicrograph().getParticles().isEmpty())
            return null;
        return getMicrograph().getParticles().get(getMicrograph().getParticles().size() - 1);
    }



    public TiltedParticle getActiveTiltedParticle()
    {
        if (hasActiveParticle()) {
            return ((UntiltedParticle) active).getTiltedParticle();
        }else {
            return null;
        }
    }



    @Override
    public TiltPairPicker getParticlePicker()
    {
        return (TiltPairPicker)picker;
    }

    // Methods

    public void zoomIn(int sx, int sy)
    {
        super.zoomIn(sx, sy);

        if(getParticlePicker().getZoom() != getFrame().getTiltedCanvas().getMagnification())
            getFrame().getTiltedCanvas().zoomIn(sx, sy);
    }

    public void zoomOut(int sx, int sy)
    {
        super.zoomOut(sx, sy);

        if(getParticlePicker().getZoom() != getFrame().getTiltedCanvas().getMagnification())
            getFrame().getTiltedCanvas().zoomOut(sx, sy);
    }

    private void addParticle(int x, int y)
    {
        try
        {
            if (active != null && getActive().getTiltedParticle() == null)
            {
                XmippMessageDialog.showInfo(frame, "Remember to pick tilted particle for each particle");
                return;
            }
            Particle tp = getMicrograph().getAlignedTiltedParticle(x, y);
            if (getMicrograph().getAddedCount() > UntiltedMicrograph.getAlignmentMin() && !getMicrograph().getTiltedMicrograph().fits(tp.getX(), tp.getY(), getParticlePicker().getSize()))
                throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("Tilted particle"));
            UntiltedParticle p = new UntiltedParticle(x, y, getMicrograph(), getParticlePicker());

            getMicrograph().addParticle(p);

            if (tp != null)
                getMicrograph().setAlignedTiltedParticle(p);
            refreshActive(p);
            frame.updateMicrographsModel();
            frame.setChanged(true);
        }
        catch (Exception e)
        {
            XmippDialog.showInfo(frame, e.getMessage());
        }

    }


    @Override
    protected void removeParticle(PickerParticle particleToRemove) {

        UntiltedParticle toRemove = (UntiltedParticle) particleToRemove;

        getMicrograph().removeParticle(toRemove);

        if (active != null && active.equals(toRemove))
        {
            if (!getMicrograph().getParticles().isEmpty())
                refreshActive(getMicrograph().getParticles().get(getMicrograph().getParticles().size() - 1));
            else
                refreshActive(null);
        }


        if (toRemove.isAdded())
            getMicrograph().initAligner();
        refresh();
        getFrame().getTiltedCanvas().repaint();
    }

    public void refreshActive(PickerParticle up)
    {

        if (up != null)
        {
            active = up;
            TiltedParticle tp = getActive().getTiltedParticle();
            if (tp != null)
            {
                int x = getFrame().getTiltedCanvas().getXOnImage(tp.getX());
                int y = getFrame().getTiltedCanvas().getYOnImage(tp.getY());

                if (!micrograph.fits(x, y, getParticlePicker().getSize()))
                    getFrame().getTiltedCanvas().moveTo(tp);
            }
        }
        else
            active = null;
        repaint();
        getFrame().getTiltedCanvas().repaint();
    }

    protected void manageActive(int x, int y)
    {

        UntiltedParticle activeParticle = getActive();

        if (!activemoved)
            return;

        if (micrograph.fits(x, y, getFrame().getParticleSize()))
        {
            moveActiveParticle(x, y);

        }
        if (activeParticle.isAdded())// added particle on matrix has been moved. Matrix
        // changed and tilted particle has to be
        // recalculated
        {
            activeParticle.setAdded(false);
            getMicrograph().initAligner();
        }
        getFrame().getTiltedCanvas().repaint();
        setActiveMoved(false);
    }

    @Override
    public void customScrollSetup(MouseEvent e){

        super.customScrollSetup(e);
        getFrame().getTiltedCanvas().setupScrollFromUntilted(e.getX(),e.getY());

    }
    @Override
    protected void doCustomPaint(Graphics2D g2)
    {
        g2.setColor(getFrame().getColor());
        for (ManualParticle p : getMicrograph().getParticles())
        {
            drawShape(g2, p, false);
        }
        if (active != null)
        {
            g2.setColor(Color.red);
            drawShape(g2, getActive(), true);
        }
        if (getFrame().drawAngles())
            drawLine(Math.toRadians(getMicrograph().getUntiltedAngle()), g2);

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

            if (frame.isEraserMode(e))
            {
                erase(e);
                getFrame().getTiltedCanvas().repaint();
                return;

            }

            UntiltedParticle activeParticle = getActive();

            if (!isDragImage(e)){

                int x = super.offScreenX(e.getX());
                int y = super.offScreenY(e.getY());

				if (active != null && !activeParticle.isAdded() && activeParticle.getTiltedParticle() != null)
					getMicrograph().addParticleToAligner(activeParticle, true);

				UntiltedParticle p = getMicrograph().getParticle(x, y, getParticlePicker().getSize());

                if(p == null && activeParticle != null && activeParticle.contains(x, y))
					p = activeParticle;
				if (p != null)
				{
					if (SwingUtilities.isLeftMouseButton(e) && e.isShiftDown())
						removeParticle(p);
					else if (SwingUtilities.isLeftMouseButton(e))
						refreshActive(p);
				}
				else if (SwingUtilities.isLeftMouseButton(e))
				{
					if (micrograph.fits(x, y, getFrame().getParticleSize()))
						addParticle(x, y);
					else
						XmippMessageDialog.showInfo(frame, XmippMessage.getOutOfBoundsMsg(String
								.format("Particle centered at %s, %s with size %s", x, y, frame.getParticlePicker().getSize())));
				}
			}
		}
	}



	@Override
	public void mouseDragged(MouseEvent e)
	{
		super.mouseDragged(e);

        // Sync screen view with tilted canvas
        if (isDragImage(e)){
            getFrame().getTiltedCanvas().syncWithUntiltedCanvas(e.getX(),e.getY());
            return;
        }

        if (frame.isPickingAvailable(e)){

            if (frame.isEraserMode(e))
			{
                erase(e);
                return;
            }

            int x = super.offScreenX(e.getX());
            int y = super.offScreenY(e.getY());

            if (active != null && micrograph.fits(x, y, getFrame().getParticleSize())){
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
			int x = super.offScreenX(e.getX());
			int y = super.offScreenY(e.getY());
			manageActive(x, y);
		}

	}
}
