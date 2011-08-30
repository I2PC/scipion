package model;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import xmipp.MDLabel;
import xmipp.MetaData;
import xmipp.Program;

public class ParticlePicker {

	private static Logger logger;
	private static String outputdir = ".";
	private static String mgselfile = "micrographs.sel";
	private static String rundir = ".";
	private static int threads;
	private static boolean fastMode;
	private static boolean incore;
	private static boolean auto;
	private static int minparticles = 100;
	private static int mintraining = 70;
	private static String trainingfn = "training.txt";
	private static String trainingmaskfn = "mask.xmp";
	private static String autofeaturesvectorfn = "auto_feature_vectors";
	private static double mincost = 0;

	private List<Family> families;
	private List<Micrograph> micrographs;
	private static ParticlePicker ppicker;

	public static Logger getLogger() {
		try {
			if (logger == null) {
				FileHandler fh = new FileHandler("PPicker.log", true);
				fh.setFormatter(new SimpleFormatter());
				logger = Logger.getLogger("PPickerLogger");
				logger.addHandler(fh);
			}
			return logger;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	public static String getXmippPath() {
		return Program.getXmippPath();
	}

	public static String getXmippPath(String relpath) {
		return getXmippPath() + File.separator + relpath;
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

	public static boolean getIsAuto() {
		return auto;
	}

	public static int getMinParticles() {
		return minparticles;
	}

	public static String getOutputDir() {
		return outputdir;
	}

	public static String getMicrographPath(String rpath) {
		return rundir + File.separator + rpath;
	}

	public static void setMicrographsSelFile(String file) {
		mgselfile = file;
	}

	public static String getMicrographsSelFile() {
		return mgselfile;
	}

	public static void setOutputDir(String dir) {
		if (!new File(dir).exists())
			throw new IllegalArgumentException("Output dir " + dir
					+ " does not exist");
		outputdir = dir;
	}

	public static String getOutputPath(String file) {
		return outputdir + File.separator + file;
	}

	public static void setThreads(int threads) {
		ParticlePicker.threads = threads;

	}

	public static void setFastMode(boolean fastMode) {
		ParticlePicker.fastMode = fastMode;
	}

	public static void setIncore(boolean incore) {
		ParticlePicker.incore = incore;
	}

	public static boolean isFastMode() {
		return fastMode;
	}

	public static boolean isIncore() {
		return incore;
	}

	public static boolean getFastMode() {
		return fastMode;
	}

	public static void setIsAuto(boolean automatic) {
		ParticlePicker.auto = automatic;

	}

	public static int getThreads() {
		return threads;
	}

	public boolean hasEmptyMicrographs(Family f) {
		for (Micrograph m : micrographs)
			if (m.getFamilyData(f).isEmpty())
				return true;
		return false;
	}

	public static FamilyState nextStep(FamilyState step) {
		if (step == FamilyState.Manual)
			return FamilyState.Supervised;
		if (step == FamilyState.Supervised)
			return FamilyState.Manual;
		return null;
	}

	private void executeProgram(String name, String args) {
		System.out.println(args);
		try {
			Program.runByName(name, args);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	private ParticlePicker() {
		this.families = new ArrayList<Family>();
		this.micrographs = new ArrayList<Micrograph>();
		loadData();
	}

	public static ParticlePicker getInstance() {
		if (ppicker == null)
			ppicker = new ParticlePicker();
		return ppicker;
	}

	public List<Family> getFamilies() {
		return families;
	}

	public List<Micrograph> getMicrographs() {
		return micrographs;
	}

	public void persistFamilies() {
		long id;
		String filename = Family.getOFilename();
		try {
			MetaData md = new MetaData();
			for (Family f : families) {
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, f.getName(), id);
				md.setValueInt(MDLabel.MDL_PICKING_COLOR,
						f.getColor().getRGB(), id);
				md.setValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, f.getSize(),
						id);
				md.setValueString(MDLabel.MDL_PICKING_FAMILY_STATE, f.getStep()
						.toString(), id);
			}
			md.write(filename);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void loadFamilyData() {
		families.clear();
		String filename = Family.getOFilename();
		if (!new File(filename).exists()) {
			families.add(Family.getDefaultFamily());
			return;
		}

		Family family;
		int rgb, size;
		FamilyState step;
		String gname;
		try {
			MetaData md = new MetaData(filename);
			long[] ids = md.findObjects();
			for (long id : ids) {
				gname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				rgb = md.getValueInt(MDLabel.MDL_PICKING_COLOR, id);
				size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
				step = FamilyState.valueOf(md.getValueString(
						MDLabel.MDL_PICKING_FAMILY_STATE, id));
				family = new Family(gname, new Color(rgb), size, step);
				families.add(family);
			}
			if (families.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No families specified on %s", filename));
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public Family getFamily(String name) {
		for (Family f : getFamilies())
			if (f.getName().equalsIgnoreCase(name))
				return f;
		return null;
	}

	public boolean existsFamilyName(String name) {
		return getFamily(name) != null;
	}

	public void loadMicrographsData() {
		String xmd = Micrograph.getIFilename();
		micrographs.clear();
		Micrograph micrograph;
		String ctf = null, filename;
		try {
			MetaData md = new MetaData(xmd);
			boolean existsctf = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			long[] ids = md.findObjects();
			for (long id : ids) {

				filename = getMicrographPath(md.getValueString(
						MDLabel.MDL_IMAGE, id));
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new Micrograph(filename, ctf);
				loadMicrographData(micrograph);
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No micrographs specified on %s", xmd));
			for (Micrograph m : micrographs)
				loadAutomaticParticles(m);

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
		}

	}

	public void loadMicrographData(Micrograph micrograph) {
		try {
			int x, y;
			Particle particle;
			String filename = micrograph.getOFilename();
			String fname;
			Family family;
			FamilyState step;
			MicrographFamilyState state;
			MicrographFamilyData mfd;
			List<MicrographFamilyData> mfdatas = new ArrayList<MicrographFamilyData>();
			if (!new File(filename).exists())
				return;
			MetaData md = new MetaData("families@" + filename);
			MetaData md2;
			for (long id : md.findObjects()) {
				fname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				state = MicrographFamilyState.valueOf(md.getValueString(
						MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, id));
				family = getFamily(fname);
				mfd = new MicrographFamilyData(micrograph, family, state);
				if (containsBlock(filename, fname)) {
					md2 = new MetaData(fname + "@" + filename);
					for (long id2 : md2.findObjects()) {

						x = md2.getValueInt(MDLabel.MDL_XINT, id2);
						y = md2.getValueInt(MDLabel.MDL_YINT, id2);
						particle = new Particle(x, y, family, micrograph);
						mfd.addManualParticle(particle);
					}
				}
				mfdatas.add(mfd);
			}
			micrograph.setFamiliesState(mfdatas);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	private boolean containsBlock(String filename, String block) {
		try {
			return Arrays.asList(MetaData.getBlocksInMetaDataFile(filename))
					.contains(block);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	public void persistMicrographsData() {
		long id;
		try {
			MetaData md;
			String block = null;
			for (Micrograph m : micrographs) {
				if (!m.hasData())
					new File(m.getOFilename()).delete();
				else {
					persistMicrographFamilies(m);
					for (MicrographFamilyData mfd : m.getFamiliesData()) {
						if (!mfd.hasManualParticles())
							continue;
						md = new MetaData();
						for (Particle p : mfd.getManualParticles()) {
							id = md.addObject();
							md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
							md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
						}
						block = mfd.getFamily().getName() + "@"
								+ m.getOFilename();
						md.writeBlock(block);
					}
				}
				persistAutomaticParticles(m);
			}

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void persistAutomaticParticles(Micrograph m) {
		if (!m.hasAutomaticParticles())
			new File(m.getAutoOFilename()).delete();
		else
			for (MicrographFamilyData mfd : m.getFamiliesData())
				persistAutomaticParticles(mfd);
	}

	public void persistMicrographFamilies(Micrograph m) {
		long id;
		try {
			MetaData md = new MetaData();
			for (MicrographFamilyData mfd : m.getFamiliesData()) {
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, mfd.getFamily()
						.getName(), id);
				md.setValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE,
						mfd.getState().toString(), id);
			}
			md.writeBlock("families@" + m.getOFilename());

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void loadAutomaticParticles(Micrograph micrograph) {
		try {
			int x, y;
			double cost;
			boolean deleted = false;
			AutomaticParticle particle;
			MetaData md;
			Family f;
			String filename = micrograph.getAutoOFilename();
			if (!new File(filename).exists())
				return;
			String[] blocks = MetaData.getBlocksInMetaDataFile(filename);
			for (String blockname : blocks) {
				md = new MetaData(blockname + "@" + filename);
				long[] ids = md.findObjects();
				for (long id : ids) {

					x = md.getValueInt(MDLabel.MDL_XINT, id);
					y = md.getValueInt(MDLabel.MDL_YINT, id);
					cost = md.getValueDouble(MDLabel.MDL_COST, id);
					f = getFamily(blockname);
					if (f == null)
						throw new IllegalArgumentException("Unknown family "
								+ blockname);
					deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? false
							: true;
					particle = new AutomaticParticle(x, y, f, micrograph, cost,
							deleted);
					micrograph.addAutomaticParticle(particle);
				}
			}
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void persistAutomaticParticles(MicrographFamilyData mfd) {
		try {
			long id;
			if (mfd.getMicrograph().hasAutomaticParticles()) {
				String filename = mfd.getMicrograph().getAutoOFilename();
				MetaData md = new MetaData();
				for (AutomaticParticle p : mfd.getAutomaticParticles()) {
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
					md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					md.setValueInt(MDLabel.MDL_ENABLED, (!p.isDeleted()) ? 1
							: -1, id);
				}
				md.write(mfd.getFamily().getName() + "@" + filename);
			}

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void removeFamily(Family family) {
		if (family.particles > 0)
			throw new IllegalArgumentException(
					Constants.getAssociatedDataMsg("family"));
		if (families.size() == 1)
			throw new IllegalArgumentException(
					Constants.getIllegalDeleteMsg("family"));
		families.remove(family);
	}

	public void train(Family f) {

		String args;
		for (Micrograph micrograph : micrographs) {
			if (!micrograph.getFamilyData(f).isEmpty()) {

				args = String
						.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train %s",
								micrograph.getFilename(),// -i
								f.getSize(), // --particleSize
								f.getOutputRoot(),// --model
								micrograph.getOutputRoot(), // --outputRoot
								f.getName() + "@" + micrograph.getOFilename());// train
																				// parameter
				if (isFastMode())
					args += " --fast";
				if (isIncore())
					args += " --in_core";
				executeProgram("xmipp_micrograph_automatic_picking", args);
			}
		}
	}

	public void classify(Family f, Micrograph micrograph) {
		String args;
		args = String
				.format("-i %s --particleSize %s --model %s --outputRoot %s --mode try --thr %s",
						micrograph.getFilename(),// -i
						f.getSize(), // --particleSize
						f.getOutputRoot(),// --model
						micrograph.getOutputRoot(),// --outputRoot
						getThreads()// --thr
				);

		if (isFastMode())
			args += " --fast";
		if (isIncore())
			args += " --in_core";
		executeProgram("xmipp_micrograph_automatic_picking", args);
	}

	public void correct(Family family, Micrograph micrograph) {
		String args;
		args = String
				.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train ",
						micrograph.getFilename(),// -i
						family.getSize(), // --particleSize
						family.getOutputRoot(),// --model
						micrograph.getOutputRoot()// --outputRoot
				);
		if (micrograph.getFamilyData(family).getManualParticles().size() > 0)
			args += family.getName() + "@" + micrograph.getOFilename();
		if (isFastMode())
			args += " --fast";
		if (isIncore())
			args += " --in_core";
		executeProgram("xmipp_micrograph_automatic_picking", args);

	}

	public int getNextFreeMicrograph(Family f) {
		int count = 0;
		for (Micrograph m : micrographs) {
			if (m.getFamilyData(f).getState() == MicrographFamilyState.Available)
				return count;
			count++;
		}
		return -1;
	}

	public void resetModel(Family family) {
		try {
			new File(family.getOTrainingFilename()).delete();
			new File(family.getOTrainingMaskFilename()).delete();
			MetaData emptymd = new MetaData();
			String block;
			MicrographFamilyData mfd;
			for (Micrograph m : micrographs) {
				mfd = m.getFamilyData(family);
			}
			saveData();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
		}
	}

	public void resetFamilyData(MicrographFamilyData mfd) {
		String block;
		MetaData emptymd = new MetaData();
		try {
			//just in case of user reset
			new File(mfd.getOTrainingAutoFeaturesVectorFilename()).delete();
			if (mfd.getStep() == FamilyState.Supervised) {
				mfd.reset();// Resetting family data
				if (mfd.getMicrograph().hasAutomaticParticles()) {
					// removing automatic particles
					block = String.format("%s@%s", mfd.getFamily().getName(),
							mfd.getMicrograph().getAutoOFilename());

					emptymd.writeBlock(block);
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void loadData() {
		loadFamilyData();
		loadMicrographsData();
	}

	public void saveData() {
		ppicker.persistFamilies();
		ppicker.persistMicrographsData();

	}
}
