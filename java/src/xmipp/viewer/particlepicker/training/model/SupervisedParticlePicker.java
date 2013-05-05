package xmipp.viewer.particlepicker.training.model;

import java.io.File;
import java.util.logging.Level;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Family;
import xmipp.viewer.particlepicker.Micrograph;



public class SupervisedParticlePicker extends TrainingPicker {
	
	private static int mintraining = 15;
	private static String trainingfn = "training.txt";
	private static String trainingmaskfn = "mask.xmp";
	private static String autofeaturesvectorfn = "auto_feature_vectors";



	
	private int threads;
	private boolean fastmode;
	private boolean incore;
	
	public SupervisedParticlePicker(String selfile, String outputdir, Integer threads,
			boolean fastmode, boolean incore) {
		this(selfile, outputdir, null, threads, fastmode, incore);
		
		

	}
	
	public SupervisedParticlePicker(String selfile, String outputdir, String fname, Integer threads,
			boolean fastmode, boolean incore) {
		super(selfile, outputdir, fname, Mode.Supervised);
		this.threads = threads;
		this.fastmode = fastmode;
		this.incore = incore;
		int count = 0;
		for (TrainingMicrograph m : micrographs)
		{
			loadMicrographData(m);
			count += m.getFamilyData(family).getParticles().size();
		}
		if(count < mintraining)
			throw new IllegalArgumentException(XmippMessage.getIllegalValueMsgWithInfo("Particles provided", count, "Must enter at least " + mintraining));
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

	public String getTrainingFile(String familyname) {
		return getOutputPath(String.format("%s_%s", familyname,
				getTrainingFilenameGeneric()));
	}

	public String getTrainingMaskFile(String familyname) {
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
	
	//not used by anyone
	public void resetModel() {
		try {
			new File(getTrainingFile(family.getName())).delete();
			new File(getTrainingMaskFile(family.getName())).delete();
			MicrographFamilyData mfd;
			for (TrainingMicrograph m : micrographs) {
				mfd = m.getFamilyData(family);
				if (mfd.getStep() == Mode.Supervised)
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

//		if (mfd.getManualParticles().size() > 0)
//			args += family.getName() + "@" + getOutputPath(micrograph.getPosFile());
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
		String args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode try --thr %s --autoPercent %s", micrograph.getFile(),// -i
		//String args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode try --thr %s", micrograph.getFile(),// -i
				family.getSize(), // --particleSize
				getOutputPath(family.getName()),// --model
				getOutputPath(micrograph.getName()),// --outputRoot
				getThreads(),//
				mfd.getAutopickpercent()// --thr
		);

		if (isFastMode())
			args += " --fast";
		if (isIncore())
			args += " --in_core";
		return args;
	}
	
	public String getTrainCommandLineArgs()
	{
		if(getManualParticlesNumber(family) < mintraining)
			throw new IllegalArgumentException(String.format("You should have at least %s particles to go to %s mode",	mintraining, Mode.Supervised));
		String args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train",
				getOutputPath(family.getName()) + "_particle_avg.xmp",//this is temporarily so it works
				family.getSize(), // --particleSize
				getOutputPath(family.getName()),// --model
				getOutputPath(family.getName()) // --outputRoot
				);// train
		// parameter
//		if (isFastMode())
//			args += " --fast";
//		if (isIncore())
//			args += " --in_core";
		return args;
	}

	public String getBuildInvariantCommandLineArgs(MicrographFamilyData mfd)
	{
		int filternum = 6;
		int NPCA = 4;
		int NCORR = 2;
		Family family = mfd.getFamily();
		Micrograph micrograph = mfd.getMicrograph();
		String args = String.format("-i %s --particleSize %s --model %s --outputRoot %s --mode buildinv %s --filter_num %s --NPCA %s --NCORR %s", micrograph.getFile(),// -i
				family.getSize(), // --particleSize
				getOutputPath(family.getName()),// --model
				getOutputPath(micrograph.getName()), // --outputRoot
				family.getName() + "@" + getOutputPath(micrograph.getPosFile()), 
				filternum,
				NPCA, 
				NCORR);// train
		// parameter
//		if (isFastMode())
//			args += " --fast";
//		if (isIncore())
//			args += " --in_core";
		return args;
	}


}
