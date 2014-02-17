package xmipp.viewer.particlepicker.training.model;

import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.ParticlePicker;

public class AutomaticParticle extends ManualParticle {
	
	private boolean deleted;

	public AutomaticParticle(int x, int y, ParticlePicker picker, SupervisedParticlePickerMicrograph micrograph, double cost, boolean deleted) {
		super(x, y, picker, micrograph, cost);
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

	public void setCost(double cost)
	{
		this.cost = cost;
		
	}
	


}
