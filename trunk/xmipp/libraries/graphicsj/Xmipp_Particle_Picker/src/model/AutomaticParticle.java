package model;

public class AutomaticParticle extends Particle {
	
	private boolean deleted;
	private double cost;

	public AutomaticParticle(int x, int y, Family family, Micrograph micrograph, double cost, boolean deleted) {
		super(x, y, family, micrograph);
		this.deleted = deleted;
		this.cost = cost;
	}
	
	public void setDeleted(boolean deleted)
	{
		this.deleted = deleted;
	}
	
	public boolean isDeleted()
	{
		return deleted;
	}
	
	public double getCost()
	{
		return cost;
	}

}
