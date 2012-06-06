package xmipp.particlepicker.training.model;

import java.io.File;
import java.util.logging.Level;

import xmipp.particlepicker.Family;
import xmipp.particlepicker.Micrograph;



public class SupervisedParticlePicker extends TrainingPicker {
	
	private static int mintraining = 70;
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

	public String getCorrectCommandLineArgs(MicrographFamilyData mfd)
	{
		Family family = mfd.getFamily();
		Micrograph micrograph = mfd.getMicrograph();
		String args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train ", micrograph.getFile(),// -i
				family.getSize(), // --particleSize
				getOutputPath(family.getName()),// --model
				getOutputPath(micrograph.getName())// --outputRoot
		);

		if (mfd.getManualParticles().size() > 0)
			args += family.getName() + "@" + getOutputPath(micrograph.getPosFile());
		if (isFastMode())
			args += " --fast";
		if (isIncore())
			args += " --in_core";
		return args;
	}
	
	public String getAutopickCommandLineArgs(MicrographFamilyData mfd)
	{
		Family family = mfd.getFamily();
		Micrograph micrograph = mfd.getMicrograph();
		String args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode try --thr %s", micrograph.getFile(),// -i
				family.getSize(), // --particleSize
				getOutputPath(family.getName()),// --model
				getOutputPath(micrograph.getName()),// --outputRoot
				getThreads()// --thr
		);

		if (isFastMode())
			args += " --fast";
		if (isIncore())
			args += " --in_core";
		return args;
	}
	
	public String getTrainCommandLineArgs(MicrographFamilyData mfd)
	{
		Family family = mfd.getFamily();
		Micrograph micrograph = mfd.getMicrograph();
		String args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train %s", micrograph.getFile(),// -i
				family.getSize(), // --particleSize
				getOutputPath(family.getName()),// --model
				getOutputPath(micrograph.getName()), // --outputRoot
				family.getName() + "@" + getOutputPath(micrograph.getPosFile()));// train
		// parameter
		if (isFastMode())
			args += " --fast";
		if (isIncore())
			args += " --in_core";
		return args;
	}


}
