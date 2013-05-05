package xmipp.viewer.particlepicker.training;

import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public class ReadOnlyParticlePicker extends TrainingPicker
{

	public ReadOnlyParticlePicker(String selfile, String outputdir)
	{
		super(selfile, outputdir, Mode.ReadOnly);
		for (TrainingMicrograph m : micrographs)
			loadMicrographData(m);
	}
	
	

}
