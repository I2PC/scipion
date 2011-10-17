package particlepicker.training.model;

import java.util.ArrayList;

import particlepicker.Family;

public class ManualParticlePicker extends TrainingPicker {
	
	public ManualParticlePicker(String selfile, String outputdir, FamilyState mode) {

		super(selfile, outputdir, mode);
		loadMicrographs();

	}



}
