package model;

import java.util.ArrayList;
import java.util.List;

public class MicrographFamilyData {
	
	private List<Particle> particles;
	private Family family;
	private Step step;
	
	

	public MicrographFamilyData(Family family) {
		this.family = family;
		this.particles = new ArrayList<Particle>();
	}
	
	public MicrographFamilyData(Family family, Step step) {
		this.family = family;
		this.particles = new ArrayList<Particle>();
		this.step = step;
	}
	public Step getStep()
	{
		return step;
	}
	
	public void setStep(Step step)
	{
		this.step = step;
	}

	public List<Particle> getParticles() {
		return particles;
	}

	
	public Family getFamily() {
		return family;
	}

	public void addParticle(Particle p) {
		particles.add(p);
		family.particles ++;
		if(step == null)
			step = family.getStep();
	}

	public void removeParticle(Particle p) {
		particles.remove(p);
		family.particles --;
		if(particles.size() == 0)
			step = null;
	}

	public boolean isEmpty()
	{
		return particles.size() == 0;
	}
	
	public boolean isReadOnly()
	{
		if(step == null)
			return false;
		return !step.equals(family.getStep());
					
	}

}
