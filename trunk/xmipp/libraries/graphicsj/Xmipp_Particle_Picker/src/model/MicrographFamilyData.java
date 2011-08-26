package model;

import java.util.ArrayList;
import java.util.List;

public class MicrographFamilyData {

	private List<Particle> particles;
	private Family family;
	private Step mystep;
	private List<AutomaticParticle> autoparticles;
	private Micrograph micrograph;
	private State state;

	public MicrographFamilyData(Micrograph micrograph, Family family) {
		this.family = family;
		this.particles = new ArrayList<Particle>();
		this.autoparticles = new ArrayList<AutomaticParticle>();
		this.micrograph = micrograph;
		state = State.Available;
		mystep = Step.Available;
	}

	public MicrographFamilyData(Micrograph micrograph, Family family,
			Step step, State state) {
		this(micrograph, family);
		this.mystep = step;
		this.state = state;
	}

	public State getState() {
		return state;
	}

	public void setState(State state) {
		if(state == State.Correct)
			this.setStep(Step.Supervised);
		this.state = state;
	}

	public Micrograph getMicrograph() {
		return micrograph;
	}

	public Step getStep() {
		return mystep;
	}

	public void setStep(Step step) {
		this.mystep = step;
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
		if (mystep == Step.Available)
			mystep = family.getStep();
		if (mystep == Step.Manual)
			state = State.Manual;

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
			if (particles.size() == 0)
				state = State.Available;
		}
	}

	public boolean isEmpty() {
		return particles.size() == 0 && autoparticles.size() == 0;
	}

	public void addAutomaticParticle(AutomaticParticle p) {
		autoparticles.add(p);

	}

	public boolean isPickingAvailable() {
		if (family.getStep() == Step.Supervised) {
			if (mystep == Step.Available)
				return false;
			if (mystep == Step.Manual)
				return false;
			if (mystep == Step.Supervised)
				if (state != State.Correct)
					return false;
			return true;
		}
		if (family.getStep() == Step.Manual) {
			if (mystep == Step.Available)
				return true;
			if (mystep == Step.Manual)
				return true;
			if (mystep == Step.Supervised)
				return false;
		}
		return true;
	}

	public boolean isActionAvailable() {
		
		if (family.getStep() == Step.Manual)
			return false;
		if (family.getStep() == Step.Supervised) {
			if (mystep == Step.Available)
				return (state == State.Available);
			if (mystep == Step.Manual)
				return false;
			if (mystep == Step.Supervised) {
				if (state == State.Manual || state == State.ReadOnly)
					return false;
			}
		}
		return true;
	}

	public String getAction() {
		if (mystep == Step.Manual)
			return null;
		if (mystep == Step.Available)
			return State.Autopick.toString();
		return state.toString();
	}

}
