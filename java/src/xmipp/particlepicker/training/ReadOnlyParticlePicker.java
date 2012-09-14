package xmipp.particlepicker.training;

import xmipp.particlepicker.training.model.FamilyState;
import xmipp.particlepicker.training.model.TrainingPicker;

public class ReadOnlyParticlePicker extends TrainingPicker
{

	public ReadOnlyParticlePicker(String selfile, String outputdir)
	{
		super(selfile, outputdir, FamilyState.ReadOnly);
		loadMicrographs();
	}
	
	

}
