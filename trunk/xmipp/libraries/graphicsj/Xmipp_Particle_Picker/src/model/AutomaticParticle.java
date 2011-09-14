package model;

public class AutomaticParticle extends Particle {
	
	private boolean deleted;

	public AutomaticParticle(int x, int y, Family family, Micrograph micrograph, double cost, boolean deleted) {
		super(x, y, family, micrograph, cost);
		if(cost> 1)
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("cost", cost));
		this.deleted = deleted;
	}
	
	public void setDeleted(boolean deleted)
	{
		this.deleted = deleted;
	}
	
	public boolean isDeleted()
	{
		return deleted;
	}
	


}
