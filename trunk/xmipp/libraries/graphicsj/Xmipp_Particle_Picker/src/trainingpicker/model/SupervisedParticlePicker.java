package trainingpicker.model;

import java.io.File;
import java.util.logging.Level;



public class SupervisedParticlePicker extends TrainingPicker {
	
	private static int mintraining = 10;
	private static String trainingfn = "training.txt";
	private static String trainingmaskfn = "mask.xmp";
	private static String autofeaturesvectorfn = "auto_feature_vectors";



	
	private int threads;
	private boolean fastmode;
	private boolean incore;
	
	public SupervisedParticlePicker(String selfile, String outputdir, Integer threads,
			boolean fastmode, boolean incore) {
		super(selfile, outputdir, FamilyState.Supervised);
		this.threads = threads;
		this.fastmode = fastmode;
		this.incore = incore;
		loadMicrographs();

	}
	
	public boolean isFastMode() {
		return fastmode;
	}

	public boolean isIncore() {
		return incore;
	}

	public int getThreads() {
		return threads;
	}
	
	
	public String getTrainingAutoFeaturesVectorFile(MicrographFamilyData mfd) {
		return getOutputPath(String.format("%s_%s_%s.txt", mfd.getMicrograph()
				.getName(), autofeaturesvectorfn, mfd.getFamily().getName()));
	}

	public String getOTrainingFilename(String familyname) {
		return getOutputPath(String.format("%s_%s", familyname,
				getTrainingFilenameGeneric()));
	}

	public String getOTrainingMaskFilename(String familyname) {
		return getOutputPath(String.format("%s_%s", familyname,
				getTrainingMaskFilenameGeneric()));
	}
	
	public static String getTrainingFilenameGeneric() {
		return trainingfn;
	}

	public static String getTrainingMaskFilenameGeneric() {
		return trainingmaskfn;
	}

	public static String getTrainingAutoFeatureVectorsFilenameGeneric() {
		return autofeaturesvectorfn;
	}

	public static int getMinForTraining() {
		return mintraining;
	}
	
	public void resetModel(Family family) {
		try {
			new File(getOTrainingFilename(family.getName())).delete();
			new File(getOTrainingMaskFilename(family.getName())).delete();
			MicrographFamilyData mfd;
			for (TrainingMicrograph m : micrographs) {
				mfd = m.getFamilyData(family);
				if (mfd.getStep() == FamilyState.Supervised)
					resetFamilyData(mfd);
			}
			saveData();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}
	


}
