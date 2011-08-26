package model;

import java.util.ArrayList;
import java.util.List;

public class MicrographFamilyData {
	
	private List<Particle> particles;
	private Family family;
	private Step step;
	private List<AutomaticParticle> autoparticles;
	private Micrograph micrograph;
	private State state;
	

	public MicrographFamilyData(Micrograph micrograph, Family family) {
		this.family = family;
		this.particles = new ArrayList<Particle>();
		this.autoparticles = new ArrayList<AutomaticParticle>();
		this.micrograph = micrograph;
		state = State.Available;
		step = Step.Available;
	}
	
	public MicrographFamilyData(Micrograph micrograph, Family family, Step step, State state) {
		this(micrograph, family);
		this.step = step;
		this.state = state;
	}
	
	public State getState() {
		return state;
	}
	public void setState(State state) {
		this.state = state;
	}
	
	public Micrograph getMicrograph()
	{
		return micrograph;
	}
	
	public Step getStep()
	{
		return step;
	}
	
	public void setStep(Step step)
	{
		this.step = step;
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
		family.particles ++;
		if(step == Step.Available)
			step = family.getStep();
		if(step == Step.Manual)
			state = State.Manual;
		
	}

	public void removeParticle(Particle p) {
		if(p == null)
			throw new IllegalArgumentException(Constants.getEmptyFieldMsg("particle"));
		if(p instanceof AutomaticParticle)
			((AutomaticParticle) p).setDeleted(true);
		else
		{
			particles.remove(p);
			family.particles --;
			if(particles.size() == 0)
				step = Step.Available;
		}
	}
	
	

	public boolean isEmpty()
	{
		return particles.size() == 0 && autoparticles.size() == 0;
	}
	

	public void addAutomaticParticle(AutomaticParticle p) {
		autoparticles.add(p);
		
	}

	public boolean isPickingAvailable() {
		if(step == Step.Available)
			return true;
		if(step == Step.Manual)
			if(!(state == State.Manual || state == State.Available))
				return false;
		if(step == Step.Supervised)
			if(state != State.Correct)
			return false;
		return true;
	}

	public boolean isActionAvailable()
	{
		if(step == Step.Available)
			return (state == State.Available && family.getStep() != Step.Manual);
		if(step == Step.Manual)
			return false;
		if(step == Step.Supervised)
		{
			if(state == State.Manual || state == State.ReadOnly)
				return false;
		}
		return true;
	}
	
	public String getAction()
	{
		if(step == Step.Manual)
			return null;
		if(step == Step.Available)
			return State.Autopic.toString();
		return state.toString();
	}

}
