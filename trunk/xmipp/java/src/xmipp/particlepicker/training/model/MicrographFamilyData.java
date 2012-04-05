package xmipp.particlepicker.training.model;

import java.util.ArrayList;
import java.util.List;

import xmipp.particlepicker.Family;
import xmipp.utils.XmippMessage;

public class MicrographFamilyData
{

	private List<TrainingParticle> manualparticles;
	private Family family;
	private List<AutomaticParticle> autoparticles;
	private TrainingMicrograph micrograph;
	private MicrographFamilyState state;

	public MicrographFamilyData(TrainingMicrograph micrograph, Family family)
	{
		if (family == null)
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("family"));
		this.family = family;
		this.manualparticles = new ArrayList<TrainingParticle>();
		this.autoparticles = new ArrayList<AutomaticParticle>();
		this.micrograph = micrograph;
		setState(MicrographFamilyState.Available);
	}

	public MicrographFamilyData(TrainingMicrograph micrograph, Family family, MicrographFamilyState state)
	{
		this(micrograph, family);
		setState(state);
	}

	public MicrographFamilyState getState()
	{
		return state;
	}

	public void setState(MicrographFamilyState state)
	{
		if (state == MicrographFamilyState.Available && !(manualparticles.isEmpty() || autoparticles.isEmpty()))
			throw new IllegalArgumentException("Micrograph has data. Can not be " + MicrographFamilyState.Available);
		this.state = state;
	}

	public TrainingMicrograph getMicrograph()
	{
		return micrograph;
	}

	public List<TrainingParticle> getManualParticles()
	{
		return manualparticles;
	}

	public List<AutomaticParticle> getAutomaticParticles()
	{
		return autoparticles;
	}

	public Family getFamily()
	{
		return family;
	}

	public void addManualParticle(TrainingParticle p)
	{

		manualparticles.add(p);
		if (state == MicrographFamilyState.Available)//to put micrograph family data on new state, done only for first particle
		{
			if (family.getStep() == FamilyState.Manual)
				state = MicrographFamilyState.Manual;
			else if (family.getStep() == FamilyState.Supervised && state == MicrographFamilyState.Autopick)//THIS STATE IS NEVER PERSISTED
				state = MicrographFamilyState.Correct;
			else if (family.getStep() == FamilyState.Review)
				state = MicrographFamilyState.Review;
			else
				throw new IllegalArgumentException(String.format("Micrograph %s could not update its state to %s and can't keep previous state %s and have particles", micrograph.getName(), state, MicrographFamilyState.Available));
		}
		
	}

	public void removeParticle(TrainingParticle p, TrainingPicker ppicker)
	{
		if (p == null)
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("particle"));
		if (p instanceof AutomaticParticle)
		{
			if (ppicker.getMode() != FamilyState.Review)
				((AutomaticParticle) p).setDeleted(true);
			else
				autoparticles.remove(p);
		}
		else
		{
			manualparticles.remove(p);
			if (manualparticles.size() == 0 && autoparticles.size() - getAutomaticParticlesDeleted() == 0)
				state = MicrographFamilyState.Available;
		}
	}

	public boolean hasManualParticles()
	{
		return !manualparticles.isEmpty();
	}

	public boolean hasAutomaticParticles()
	{
		return autoparticles.size() != 0;
	}

	public boolean isEmpty()
	{
		return manualparticles.isEmpty() && autoparticles.isEmpty();
	}

	public int getAutomaticParticles(double threshold)
	{
		if (autoparticles.isEmpty())
			return 0;
		return autoparticles.size() - getAutomaticParticlesDeleted(threshold);
	}

	public void addAutomaticParticle(AutomaticParticle p)
	{
		if (state == MicrographFamilyState.Available)
			throw new IllegalArgumentException(String.format("Invalid state %s on micrograph %s and family %s for adding automatic particles", state, micrograph.getName(), family.getName()));
		autoparticles.add(p);

	}

	public boolean isPickingAvailable()
	{
		if (family.getStep() == FamilyState.Supervised)
		{
			if (state == MicrographFamilyState.Available)
				return false;
			if (state == MicrographFamilyState.Manual)
				return false;
			if (state != MicrographFamilyState.Correct)
				return false;
			return true;
		}
		if (family.getStep() == FamilyState.Manual)
		{
			if (state == MicrographFamilyState.Available)
				return true;
			if (state == MicrographFamilyState.Manual)
				return true;
			return false;
		}
		if (family.getStep() == FamilyState.Review)
			return true;
		return false;
	}

	public boolean isActionAvailable(double threshold)
	{

		if (family.getStep() != FamilyState.Supervised)
			return false;
		if (family.getStep() == FamilyState.Supervised)
		{
			if (state == MicrographFamilyState.Available)
				return true;
			if (state == MicrographFamilyState.Manual)
				return false;
			if (state == MicrographFamilyState.ReadOnly)
				return false;
			if (state == MicrographFamilyState.Correct)
				//return !manualparticles.isEmpty() || getAutomaticParticlesDeleted(threshold) > 0;
				return true;
			return true;
		}
		return false;
	}

	public String getAction()
	{
		if (state == MicrographFamilyState.Manual)
			return null;
		if (state == MicrographFamilyState.Available)
			return MicrographFamilyState.Autopick.toString();
		return state.toString();
	}

	public void reset()
	{
		autoparticles.clear();
		manualparticles.clear();
		setState(MicrographFamilyState.Available);
	}

	public FamilyState getStep()
	{
		if (state == MicrographFamilyState.Manual)
			return FamilyState.Manual;
		if (state == MicrographFamilyState.Available)
			return FamilyState.Available;
		return FamilyState.Supervised;
	}

	public int getAutomaticParticlesDeleted()
	{
		int count = 0;
		for (AutomaticParticle p : autoparticles)
			if (p.isDeleted())
				count++;
		return count;
	}

	public int getAutomaticParticlesCount(double threshold)
	{
		return autoparticles.size() - getAutomaticParticlesDeleted(threshold);
	}

	public int getAutomaticParticlesDeleted(double threshold)
	{
		int count = 0;
		for (AutomaticParticle p : autoparticles)
			if (p.isDeleted() || p.getCost() < threshold)
				count++;
		return count;
	}

	public void deleteBelowThreshold(double threshold)
	{
		for (AutomaticParticle p : autoparticles)
			if (p.getCost() < threshold)
				p.setDeleted(true);
	}

	public List<TrainingParticle> getParticles()
	{
		ArrayList<TrainingParticle> result = new ArrayList<TrainingParticle>();
		result.addAll(manualparticles);
		result.addAll(autoparticles);
		return result;
	}
	
	public List<TrainingParticle> getAvailableParticles(double threshold)
	{
		ArrayList<TrainingParticle> result = new ArrayList<TrainingParticle>();
		result.addAll(manualparticles);
		for(AutomaticParticle ap: autoparticles)
			if(!ap.isDeleted() && ap.getCost() >= threshold)
				result.add(ap);
		return result;
	}
	
	public TrainingParticle getLastAvailableParticle()
	{
		AutomaticParticle ap;
		for(int i = autoparticles.size() - 1; i >= 0; i --)
		{
			ap = autoparticles.get(i);
			if(!ap.isDeleted())
				return ap;
		}
		if(!manualparticles.isEmpty())
			return manualparticles.get(manualparticles.size() - 1);
		return null;
		
	}
	
	public String toString()
	{
		return String.format("Micrograph: %s Family: %s State: %s", micrograph.getName(), family.getName(), state);
	}

}
