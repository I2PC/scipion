package xmipp.viewer.particlepicker.training.model;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import xmipp.jni.Particle;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.PickerParticle;
import java.awt.Rectangle;

public class SupervisedPickerMicrograph extends Micrograph
{

	private boolean autopicking = false;	
	private List<ManualParticle> manualparticles;
	private List<AutomaticParticle> autoparticles;
	private MicrographState state;
	private int autopickpercent;
	private double threshold = 0.0;
    private Rectangle rectangle;

	public double getThreshold() {
		return threshold;
	}



	public void setThreshold(double threshold) {
		this.threshold = threshold;
	}



	public SupervisedPickerMicrograph(String file, String psd, String ctf)
	{
		super(file, psd, ctf);
                
		this.manualparticles = new ArrayList<ManualParticle>();
		this.autoparticles = new ArrayList<AutomaticParticle>();
		state = MicrographState.Available;
	}



	public boolean isAutopicking()
	{
		return autopicking;
	}

	public void setAutopicking(boolean autopicking)
	{
		this.autopicking = autopicking;
	}

	public static String getName(String file)
	{
		String[] tokens = file.split(File.separator);
		if (tokens.length < 2)
			throw new IllegalArgumentException("Name for micrograph" + "is taken from parent dir, invalid path " + file);
		return tokens[tokens.length - 2];
	}



	public ManualParticle getParticle(int x, int y)
	{
		for (ManualParticle p : getManualParticles())
			if (p.contains(x, y))
				return p;
		return null;
	}

	public AutomaticParticle getAutomaticParticle(int x, int y, double threshold)
	{
		for (AutomaticParticle p : getAutomaticParticles())
			if (!p.isDeleted() && p.getCost() >= threshold && p.contains(x, y))
				return p;
		return null;
	}

	
	public boolean hasManualParticles()
	{
		if (getManualParticles().size() > 0)
			return true;
		return false;

	}

	public void removeParticles(int x, int y, SupervisedParticlePicker ppicker)
	{
		List<ManualParticle> particles = new ArrayList<ManualParticle>();

		for (ManualParticle p : getManualParticles())
			if (p.contains(x, y))
				particles.add(p);
		for (ManualParticle p : particles)
			removeParticle(p, ppicker);
		particles.clear();
		for (AutomaticParticle p : getAutomaticParticles())
			if (p.contains(x, y))
				particles.add(p);

		for (ManualParticle p : particles)
			removeParticle(p, ppicker);

	}

	public List<ManualParticle> getManualParticles()
	{
		return manualparticles;
	}

	public List<AutomaticParticle> getAutomaticParticles()
	{
		return autoparticles;
	}

	public void addManualParticle(ManualParticle p, SupervisedParticlePicker ppicker)
	{
		manualparticles.add(p);

		if (state == MicrographState.Available || state == MicrographState.Auto)
		{
			if (ppicker.getMode() == Mode.Manual)
				state = MicrographState.Manual;
			else if (ppicker.getMode() == Mode.Review)
				state = MicrographState.Review;
			else
				throw new IllegalArgumentException(
						String.format("Micrograph %s could not update its state to %s and can't keep previous state %s and have particles", getName(), state, MicrographState.Available));
		}
		
		
	}

	public void removeParticle(PickerParticle p, SupervisedParticlePicker ppicker)
	{
		if (p == null)
			throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("particle"));
		if (p instanceof AutomaticParticle)
		{
			if (ppicker.getMode() != Mode.Review)
				((AutomaticParticle) p).setDeleted(true);
			else
				autoparticles.remove(p);
			// deleting could be the first thing to do after autopick, so I have
			// to mark micrograph on this choice too
			if (state == MicrographState.Auto && ppicker.getMode() == Mode.Review)
				state = MicrographState.Review;
		}
		else
			manualparticles.remove(p);
		

	}

	public boolean hasAutomaticParticles()
	{
		return autoparticles.size() != 0;
	}

	public boolean isEmpty()
	{
		return manualparticles.isEmpty() && autoparticles.isEmpty();
	}

	public int getAutomaticParticlesNumber(double threshold)
	{
		if (autoparticles.isEmpty())
			return 0;
		return autoparticles.size() - getAutomaticParticlesDeleted(threshold);
	}
        
        public int getAutomaticParticlesNumber()
	{
		if (autoparticles.isEmpty())
			return 0;
		return autoparticles.size() - getAutomaticParticlesDeleted();
	}

	public synchronized void addAutomaticParticle(AutomaticParticle p)
	{
		autoparticles.add(p);

	}

	

	public MicrographState getState()
	{
		return state;
	}

	public void setState(MicrographState state)
	{
		if (state == MicrographState.Available && !(manualparticles.isEmpty() || autoparticles.isEmpty()))
			throw new IllegalArgumentException("Micrograph has data. Can not be " + MicrographState.Available);
		this.state = state;
	}

	public Mode getStep()
	{
		if (state == MicrographState.Manual)
			return Mode.Manual;
		if (state == MicrographState.Available)
			return Mode.Available;
		return Mode.Supervised;
	}

	public int getAutomaticParticlesDeleted()
	{
		int count = 0;
		for (AutomaticParticle p : autoparticles)
			if (p.isDeleted())
				count++;
		return count;
	}

	

	public synchronized int getAutomaticParticlesDeleted(double threshold)
	{
		int count = 0;
		for (AutomaticParticle p : autoparticles)
			if (p.isDeleted() || p.getCost() < threshold)
				count++;
		return count;
	}


	public List<ManualParticle> getParticles()
	{
		ArrayList<ManualParticle> result = new ArrayList<ManualParticle>();
		result.addAll(manualparticles);
		result.addAll(autoparticles);
		return result;
	}

	public List<ManualParticle> getAvailableParticles(double threshold)
	{
		ArrayList<ManualParticle> result = new ArrayList<ManualParticle>();
		result.addAll(manualparticles);
		for (AutomaticParticle ap : autoparticles)
			if (!ap.isDeleted() && ap.getCost() >= threshold)
				result.add(ap);
		return result;
	}

	public ManualParticle getLastAvailableParticle(double threshold)
	{
		AutomaticParticle ap;
		for (int i = autoparticles.size() - 1; i >= 0; i--)
		{
			ap = autoparticles.get(i);
			if (!ap.isDeleted() && ap.getCost() >= threshold)
				return ap;
		}
		if (!manualparticles.isEmpty())
			return manualparticles.get(manualparticles.size() - 1);
		return null;

	}

	

	@Override
	public boolean hasData()
	{
		return !isEmpty();
	}

    @Override
    public List<? extends PickerParticle> getParticleList() {
        return getParticles();
    }

    public int getAutopickpercent()
	{
		return autopickpercent;
	}

	public void setAutopickpercent(int autopickpercent)
	{
		this.autopickpercent = autopickpercent;
	}
	
	public void initParticlesRectangle(SupervisedParticlePicker picker)
	{
		double x1 = Double.POSITIVE_INFINITY, y1 = Double.POSITIVE_INFINITY, x2 = Double.NEGATIVE_INFINITY, y2 = Double.NEGATIVE_INFINITY;
		List<ManualParticle> particles = getParticles();
		if(particles.isEmpty())
			return;
		for(Particle p: getParticles())
		{
			if(p.getX() < x1)
				x1 = p.getX();
			if(p.getY() < y1)
				y1 = p.getY();
			if(p.getX() > x2)
				x2 = p.getX();
			if(p.getY() > y2)
				y2 = p.getY();
		}
		
		int radius = picker.getRadius();
		radius = radius + radius/5;
		int x = Math.max(0, (int)(x1 - radius));
		int y = Math.max(0, (int)(y1 - radius));
		int width  = Math.min((int)(x2 - x + radius), this.getWidth() - 1);
		int height = Math.min((int)(y2 - y + radius), this.getHeigth() - 1);
                rectangle = new Rectangle(x, y, width, height);
		return;
		
	}
        
        public void resetParticlesRectangle()
        {
            rectangle = null;
        }
        
        public boolean isValid(SupervisedParticlePicker picker)
        {
                int count = getParticles().size();
                if(count <= 3)
                    return false;
                return true;
        
        }
	
	public void deleteBelowThreshold() {
		// TODO Auto-generated method stub
		for (AutomaticParticle p : autoparticles) 
			if(!p.isDeleted())
				p.setDeleted(p.getCost() < getThreshold());
	}

        public Rectangle getRectangle() {
            return rectangle;
        }
        
        

}
