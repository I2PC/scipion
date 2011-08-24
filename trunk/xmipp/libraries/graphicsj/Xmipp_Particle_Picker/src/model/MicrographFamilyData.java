package model;

import java.util.ArrayList;
import java.util.List;

public class MicrographFamilyData {
	
	private List<Particle> particles;
	private Family family;

	public MicrographFamilyData(Family family) {
		this.family = family;
		this.particles = new ArrayList<Particle>();
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
	}

	public void removeParticle(Particle p) {
		particles.remove(p);
		family.particles --;
	}


	
	

}
