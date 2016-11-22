package xmipp.viewer.particlepicker.tiltpair.model;

import java.util.ArrayList;
import java.util.List;
import xmipp.jni.Particle;
import xmipp.jni.TiltPairAligner;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.PickerParticle;

public class UntiltedMicrograph extends Micrograph
{

	private static int alignmentmin = 4;
	private TiltedMicrograph tiltedmicrograph;
	private List<UntiltedParticle> particles;
	private TiltPairAligner tpa;
	private int added;
	private double[] angles;

	public UntiltedMicrograph(String file, TiltedMicrograph tiltedmicrograph)
	{
		super(file, getName(file, 1));
		this.tiltedmicrograph = tiltedmicrograph;
		particles = new ArrayList<UntiltedParticle>();
		tpa = new TiltPairAligner();

	}

	public int getAddedCount()
	{
		return added;
	}

	public TiltedMicrograph getTiltedMicrograph()
	{
		return tiltedmicrograph;
	}

	public UntiltedParticle getParticle(int x, int y, int size)
	{
		for (UntiltedParticle p : particles)
			if (p.contains(x, y, size))
				return p;
		return null;

	}

	@Override
	public boolean hasData()
	{
		return !particles.isEmpty();
	}

    @Override
    public List<? extends PickerParticle> getParticleList() {
        return getParticles();
    }


    public void reset(TiltPairPicker picker)
	{
		angles = null;
		for(UntiltedParticle p: particles)
			p.setTiltedParticle(null);
		getParticles().clear();
		getTiltedMicrograph().getParticles().clear();
		//JMRT: I think is very very dangerous to delete the existing pos files.
		//new File(picker.getOutputPath(getPosFile())).delete();
		//new File(picker.getOutputPath(getTiltedMicrograph().getPosFile())).delete();
		initAligner();
	}

	public void removeParticle(UntiltedParticle p)
	{
		particles.remove(p);
		tiltedmicrograph.removeParticle(p.getTiltedParticle());
	}

	public void addParticle(UntiltedParticle p)
	{
		particles.add(p);
	}

	public List<UntiltedParticle> getParticles()
	{
		return particles;
	}

	public void addParticleToAligner(UntiltedParticle up, boolean recalculateAngles)
	{
		if (up.getTiltedParticle() == null)
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("TiltedParticle"));
		tpa.addParticleToAligner(up.getX(), up.getY(), up.getTiltedParticle().getX(), up.getTiltedParticle().getY());
		up.setAdded(true);
		added++;
		if (anglesAvailable() && recalculateAngles)
			angles = tpa.computeAngles();
	}

	public void initAligner()
	{
		tpa.clear();
		added = 0;
		for (UntiltedParticle p : particles)
			if (p.getTiltedParticle() != null)
				addParticleToAligner(p, false);
		if (anglesAvailable())
			angles = tpa.computeAngles();
	}

	public boolean anglesAvailable()
	{
		return !(added < 10);
	}

	//	public double[] getAngles()
	//	{
	//		return angles;
	//	}

	public double getUntiltedAngle()
	{
		if (angles == null)
			return 90;
		return angles[0];
	}

	public double getTiltedAngle()
	{
		if (angles == null)
			return 90;
		return angles[1];
	}

	public double getTiltAngle()
	{
		if (angles == null)
			return 0;
		return angles[2];
	}

	public Particle getAlignedTiltedParticle(int x, int y)
	{
		if (getAddedCount() < getAlignmentMin())
			return null;
		Particle p = tpa.getTiltedParticle(x, y);
		return p;
	}
	
	public Particle getAlignedUntiltedParticle(int x, int y)
	{
		if (getAddedCount() < getAlignmentMin())
			return null;
		Particle p = tpa.getUntiltedParticle(x, y);
		return p;
	}

	public void setAlignedTiltedParticle(UntiltedParticle up)
	{
        tiltedmicrograph.removeParticle(up.getTiltedParticle());
		Particle p = getAlignedTiltedParticle(up.getX(), up.getY());
        if(!tiltedmicrograph.fits(p.getX(), p.getY(), up.getParticlePicker().getSize()))
                p = null;
		if (p != null)
		{
			TiltedParticle tp = new TiltedParticle(p.getX(), p.getY(), up);
			getTiltedMicrograph().addParticle(tp);
			up.setTiltedParticle(tp);
		}
	}

	public void removeParticles(int x, int y)
	{
		List<UntiltedParticle> particles = new ArrayList<UntiltedParticle>();
		for (UntiltedParticle p : getParticles())
			if (p.contains(x, y, p.getParticlePicker().getSize()))
				particles.add(p);
		for (UntiltedParticle p : particles)
			removeParticle(p);
	}

	public static int getAlignmentMin()
	{
		return alignmentmin;
	}

}
