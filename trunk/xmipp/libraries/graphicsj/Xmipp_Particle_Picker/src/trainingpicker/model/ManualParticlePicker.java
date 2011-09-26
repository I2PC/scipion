package trainingpicker.model;

import java.util.ArrayList;

public class ManualParticlePicker extends TrainingPicker {
	
	public ManualParticlePicker(String selfile, String outputdir, FamilyState mode) {

		super(selfile, outputdir, mode);
		loadMicrographs();

	}

}
