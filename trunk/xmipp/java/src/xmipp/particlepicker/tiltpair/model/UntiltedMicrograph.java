package xmipp.particlepicker.tiltpair.model;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;

import xmipp.particlepicker.Micrograph;
import xmipp.utils.XmippMessage;
import xmipp.jni.Particle;
import xmipp.jni.TiltPairAligner;

public class UntiltedMicrograph extends Micrograph
{

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
	public void reset()
	{

		getParticles().clear();
		getTiltedMicrograph().reset();
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

	public void addParticleToAligner(UntiltedParticle up)
	{
		if (up.getTiltedParticle() == null)
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("TiltedParticle"));
		tpa.addParticleToAligner(up.getX(), up.getY(), up.getTiltedParticle().getX(), up.getTiltedParticle().getY());
		up.setAdded(true);
		added++;
		if (anglesAvailable())
		{
			angles = tpa.computeAngles();
			System.out.printf("untilted: %.2f tilted: %.2f deviation: %.2f\n", angles[0], angles[1], angles[2]);
		}
	}

	public void initAligner()
	{

		tpa.clear();
		added = 0;
		for (UntiltedParticle p : particles)
			if (p.getTiltedParticle() != null)
				addParticleToAligner(p);

	}

	public boolean anglesAvailable()
	{
		return !(added < 4);
	}

//	public double[] getAngles()
//	{
//		return angles;
//	}
	
	public double getUntiltedAngle()
	{
		if(angles == null)
			return 0;
		return angles[0];
	}
	
	public double getTiltedAngle()
	{
		if(angles == null)
			return 0;
		return angles[1];
	}
	
	public double getTiltAngle()
	{
		if(angles == null)
			return 0;
		return angles[2];
	}

	public Particle getAlignerTiltedParticle(int x, int y)
	{
		Particle p = tpa.getTiltedParticle(x, y);
		return p;
	}

	public void setAlignerTiltedParticle(UntiltedParticle up)
	{
		Particle p = getAlignerTiltedParticle(up.getX(), up.getY());

		TiltedParticle tp = new TiltedParticle(p.getX(), p.getY(), up);
		getTiltedMicrograph().addParticle(tp);
		up.setTiltedParticle(tp);
	}

}
