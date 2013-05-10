package xmipp.viewer.particlepicker.training.model;

import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.ParticlePicker;

public class AutomaticParticle extends TrainingParticle {
	
	private boolean deleted;

	public AutomaticParticle(int x, int y, ParticlePicker picker, TrainingMicrograph micrograph, double cost, boolean deleted) {
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
	


}
