package xmipp.particlepicker.training.model;

import xmipp.particlepicker.Family;
import xmipp.utils.XmippMessage;

public class AutomaticParticle extends TrainingParticle {
	
	private boolean deleted;

	public AutomaticParticle(int x, int y, Family family, TrainingMicrograph micrograph, double cost, boolean deleted) {
		super(x, y, family, micrograph, cost);
		if(cost> 1)
			throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("cost", cost));
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
