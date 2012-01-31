package xmipp.particlepicker.tiltpair.model;

import java.util.ArrayList;
import java.util.List;

import xmipp.particlepicker.Constants;
import xmipp.particlepicker.Micrograph;
import xmipp.jni.Particle;
import xmipp.jni.TiltPairAligner;

public class UntiltedMicrograph extends Micrograph {

	private TiltedMicrograph tiltedmicrograph;
	private List<UntiltedParticle> particles;
	private TiltPairAligner tpa;
	private int added;

	public UntiltedMicrograph(String file, TiltedMicrograph tiltedmicrograph) {
		super(file, getName(file, 1));
		this.tiltedmicrograph = tiltedmicrograph;
		particles = new ArrayList<UntiltedParticle>();
		tpa = new TiltPairAligner();
		
	}

	public int getAddedCount()
	{
		return added;
	}

	

	public TiltedMicrograph getTiltedMicrograph() {
		return tiltedmicrograph;
	}

	public UntiltedParticle getParticle(int x, int y, int size) {
		for (UntiltedParticle p : particles)
			if (p.contains(x, y, size))
				return p;
		return null;

	}

	@Override
	public boolean hasData() {
		return !particles.isEmpty();
	}

	@Override
	public void reset() {

		getParticles().clear();
		getTiltedMicrograph().reset();
		initAligner();
	}

	public void removeParticle(UntiltedParticle p) {
		particles.remove(p);
		tiltedmicrograph.removeParticle(p.getTiltedParticle());
	}

	public void addParticle(UntiltedParticle p) {
		particles.add(p);
	}
	
	public List<UntiltedParticle> getParticles()
	{
		return particles;
	}
	

	
	public void addParticleToAligner(UntiltedParticle up)
	{
		if(up.getTiltedParticle() == null)
			throw new IllegalArgumentException(Constants.getEmptyFieldMsg("TiltedParticle"));
		System.out.printf("Adding Particle x1 = %s y1 = %s x2 = %s y2 = %s\n", up.getX(), up.getY(), up.getTiltedParticle().getX(), up.getTiltedParticle().getY());
		tpa.addParticleToAligner(up.getX(), up.getY(), up.getTiltedParticle().getX(), up.getTiltedParticle().getY());
		up.setAdded(true);
		added ++;
		
	}
	
	
	public void initAligner()
	{
		
		tpa.clear();
		System.out.println("Initializing aligner. Already cleared ...");
		added = 0;
		for(UntiltedParticle p: particles)
			if(p.getTiltedParticle() != null)
				addParticleToAligner(p);
		System.out.println("Initialized aligner ...");
	}
	
	public int[] getAngles()
	{
		
		if(added < 4)
			throw new IllegalArgumentException("Not enough particles yet");
		System.out.println("Calling getAngles...");
		int[] angles = tpa.computeAngles();
		for(int i = 0; i < angles.length; i ++)
			System.out.println(angles[i]);
		return tpa.computeAngles();
	}
	
	public void setAlignerTiltedParticle(UntiltedParticle up)
	{
		System.out.println("Calling getTiltedParticle...");
		Particle p = tpa.getTiltedParticle(up.getX(), up.getY()); 
		System.out.println(p);
		TiltedParticle tp = new TiltedParticle(p.getX(), p.getY(), up);
		getTiltedMicrograph().addParticle(tp);
		up.setTiltedParticle(tp);
	}


}
