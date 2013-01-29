package xmipp.viewer.particlepicker.training;

import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public class ReadOnlyParticlePicker extends TrainingPicker
{

	public ReadOnlyParticlePicker(String selfile, String outputdir)
	{
		super(selfile, outputdir, FamilyState.ReadOnly);
		for (TrainingMicrograph m : micrographs)
			loadMicrographData(m);
	}
	
	

}
