package xmipp.particlepicker.training.model;

import java.util.ArrayList;

import xmipp.particlepicker.Family;

public class ManualParticlePicker extends TrainingPicker {
	
	public ManualParticlePicker(String selfile, String outputdir, FamilyState mode) {

		super(selfile, outputdir, mode);
		loadMicrographs();

	}
	
	public ManualParticlePicker(String selfile, String outputdir, String fname, FamilyState mode) {

		super(selfile, outputdir, fname, mode);
		loadMicrographs();

	}




}
