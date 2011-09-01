package model;

import java.util.ArrayList;
import java.util.List;

public class MicrographFamilyData {

	private List<Particle> particles;
	private Family family;
	private List<AutomaticParticle> autoparticles;
	private Micrograph micrograph;
	private MicrographFamilyState state;

	public MicrographFamilyData(Micrograph micrograph, Family family) {
		this.family = family;
		this.particles = new ArrayList<Particle>();
		this.autoparticles = new ArrayList<AutomaticParticle>();
		this.micrograph = micrograph;
		setState(MicrographFamilyState.Available);
	}

	public MicrographFamilyData(Micrograph micrograph, Family family,
			MicrographFamilyState state) {
		this(micrograph, family);
		
		this.state = state;
	}

	public MicrographFamilyState getState() {
		return state;
	}

	public void setState(MicrographFamilyState state) {
		if(state == MicrographFamilyState.Available && 
			!(particles.isEmpty() || autoparticles.isEmpty()))
			throw new IllegalArgumentException("Micrograph has data. Can not be " + MicrographFamilyState.Available);
		this.state = state;
	}

	public Micrograph getMicrograph() {
		return micrograph;
	}



	public List<Particle> getManualParticles() {
		return particles;
	}

	public List<AutomaticParticle> getAutomaticParticles() {
		return autoparticles;
	}

	public Family getFamily() {
		return family;
	}

	public void addManualParticle(Particle p) {

		particles.add(p);
		family.particles++;
		if (state == MicrographFamilyState.Available)
		{
			if(family.getStep() == FamilyState.Manual)
				state = MicrographFamilyState.Manual;
			else if(family.getStep() == FamilyState.Supervised && state == MicrographFamilyState.Autopick)
				state = MicrographFamilyState.Correct;
			else
				throw new IllegalArgumentException(String.format("Micrograph could not update its state and can't keep previous state %s and have particles", MicrographFamilyState.Available));
		}
	}

	
	public void removeParticle(Particle p) {
		if (p == null)
			throw new IllegalArgumentException(
					Constants.getEmptyFieldMsg("particle"));
		if (p instanceof AutomaticParticle)
			((AutomaticParticle) p).setDeleted(true);
		else {
			particles.remove(p);
			family.particles--;
			if (particles.size() == 0 && autoparticles.size() - getAutomaticParticlesDeleted() == 0)
				state = MicrographFamilyState.Available;
		}
	}
	
	public boolean hasManualParticles()
	{
		return particles.size() != 0;
	}
	
	public boolean hasAutomaticParticles()
	{
		return autoparticles.size() != 0;
	}

	public boolean isEmpty() {
		return particles.size() == 0 && autoparticles.size() == 0;
	}
	
	public int getAutomaticParticles(double threshold)
	{
		return autoparticles.size() - getAutomaticParticlesDeleted(threshold);
	}

	public void addAutomaticParticle(AutomaticParticle p) {
		if(state == MicrographFamilyState.Available)
			throw new IllegalArgumentException(String.format("Invalid state %s on micrograph %s and family %s for adding automatic particles", state, micrograph.getName(), family.getName()));
		autoparticles.add(p);

	}

	public boolean isPickingAvailable() {
		if (family.getStep() == FamilyState.Supervised) {
			if (state == MicrographFamilyState.Available)
				return false;
			if (state == MicrographFamilyState.Manual)
				return false;
			if (state != MicrographFamilyState.Correct)
					return false;
			return true;
		}
		if (family.getStep() == FamilyState.Manual) {
			if (state == MicrographFamilyState.Available)
				return true;
			if (state == MicrographFamilyState.Manual)
				return true;
			return false;
		}
		return false;
	}

	public boolean isActionAvailable(double threshold) {
		
		if (family.getStep() == FamilyState.Manual)
			return false;
		if (family.getStep() == FamilyState.Supervised) {
			if (state == MicrographFamilyState.Available)
				return true;
			if (state == MicrographFamilyState.Manual)
				return false;
			if (state == MicrographFamilyState.ReadOnly)
				return false;
			if(state == MicrographFamilyState.Correct)
				return !particles.isEmpty() || getAutomaticParticlesDeleted(threshold) > 0;
			return true;
		}
		return false;
	}

	public String getAction() {
		if (state == MicrographFamilyState.Manual)
			return null;
		if (state == MicrographFamilyState.Available)
			return MicrographFamilyState.Autopick.toString();
		return state.toString();
	}
	
	public void reset()
	{
		setState(MicrographFamilyState.Available);
		family.particles -= particles.size();
		autoparticles.clear();
		particles.clear();
	}
	
	public FamilyState getStep()
	{
		if(state == MicrographFamilyState.Manual)
			return FamilyState.Manual;
		if(state == MicrographFamilyState.Available)
			return FamilyState.Available;
		return FamilyState.Supervised;
	}
	
	public int getAutomaticParticlesDeleted()
	{
		int count = 0;
		for(AutomaticParticle p: autoparticles)
			if(p.isDeleted())
				count ++;
		return count;
	}
	
	public int getAutomaticParticlesCount(double threshold)
	{
		return autoparticles.size() - getAutomaticParticlesDeleted(threshold);
	}
	
	public int getAutomaticParticlesDeleted(double threshold)
	{
		int count = 0;
		for(AutomaticParticle p: autoparticles)
			if(p.isDeleted() || p.getCost() < threshold)
				count ++;
		return count;
	}
	
	public String getOTrainingAutoFeaturesVectorFilename()
	{
		return ParticlePicker.getOutputPath(String.format("%s_%s_%s.txt", micrograph.getName(), ParticlePicker.getTrainingAutoFeatureVectorsFilenameGeneric(), family.getName()));
	}
	
	public void deleteBelowThreshold(double threshold)
	{
		for(AutomaticParticle p: autoparticles)
			if(p.getCost() < threshold)
				p.setDeleted(true);
	}
	
}
