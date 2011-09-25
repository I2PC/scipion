package trainingpicker.model;


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
import xmipp.Particle;

public abstract class ParticlePicker {

	private String familiesfile;
	private static Logger logger;
	private String outputdir = ".";
	private static String rundir = ".";
	private String selfile = "micrographs.sel";
	
	private FamilyState mode;
	private boolean changed;
	protected List<Family> families;
	protected List<TrainingMicrograph> micrographs;

	public static FamilyState previousStep(FamilyState step) {
		if (step == FamilyState.Manual)
			return null;
		if (step == FamilyState.Supervised)
			return FamilyState.Manual;
		return null;
	}

	public ParticlePicker(String selfile, String outputdir, FamilyState mode) {

		this.selfile = selfile;
		this.outputdir = outputdir;
		this.mode = mode;
		this.families = new ArrayList<Family>();
		this.micrographs = new ArrayList<TrainingMicrograph>();
		this.familiesfile = getOutputPath("families.xmd");
		loadFamilies();

	}


	public void setChanged(boolean changed) {
		this.changed = changed;
	}

	public boolean isChanged() {
		return changed;
	}

	

	public FamilyState getMode() {
		return mode;
	}



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

	public String getOutputPath(String file) {
		return outputdir + File.separator + file;
	}

	public static String getXmippPath() {
		return Program.getXmippPath();
	}

	public static String getXmippPath(String relpath) {
		return getXmippPath() + File.separator + relpath;
	}

	public String getOutputDir() {
		return outputdir;
	}

	public static String getMicrographPath(String rpath) {
		return rundir + File.separator + rpath;
	}

	public String getMicrographsSelFile() {
		return selfile;
	}

	public boolean hasEmptyMicrographs(Family f) {
		for (TrainingMicrograph m : micrographs)
			if (m.getFamilyData(f).isEmpty())
				return true;
		return false;
	}

	public static FamilyState nextStep(FamilyState step) {
		if (step == FamilyState.Manual)
			return FamilyState.Supervised;
		if (step == FamilyState.Supervised)
			return FamilyState.Review;
		return null;
	}

	public List<Family> getFamilies() {
		return families;
	}

	public List<TrainingMicrograph> getMicrographs() {
		return micrographs;
	}

	public void persistFamilies() {
		long id;
		String filename = familiesfile;
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

	public void loadFamilies() {
		families.clear();
		String filename = familiesfile;
		if (!new File(filename).exists()) {
			families.add(Family.getDefaultFamily());
			return;
		}

		Family family;
		int rgb, size;
		FamilyState step;
		String name;
		try {
			MetaData md = new MetaData(filename);
			long[] ids = md.findObjects();
			for (long id : ids) {
				name = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				rgb = md.getValueInt(MDLabel.MDL_PICKING_COLOR, id);
				size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
				step = FamilyState.valueOf(md.getValueString(
						MDLabel.MDL_PICKING_FAMILY_STATE, id));
				if (getMode() == FamilyState.Review
						&& step != FamilyState.Review) {
					step = FamilyState.Review;
					setChanged(true);
				}
				if(step == FamilyState.Supervised && this instanceof SupervisedParticlePicker)
					if(!new File(((SupervisedParticlePicker)this).getOTrainingFilename(name)).exists())
						throw new IllegalArgumentException(String.format("Training file does not exist. Family can not be in %s mode", step));
				family = new Family(name, new Color(rgb), size, step, this);
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

	public TrainingMicrograph getMicrograph(String name) {
		for (TrainingMicrograph m : getMicrographs())
			if (m.getName().equalsIgnoreCase(name))
				return m;
		return null;
	}

	public boolean existsFamilyName(String name) {
		return getFamily(name) != null;
	}


	public void loadMicrographs() {

		micrographs.clear();
		TrainingMicrograph micrograph;
		String ctf = null, file;
		try {
			MetaData md = new MetaData(selfile);
			boolean existsctf = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			long[] ids = md.findObjects();
			for (long id : ids) {

				file = getMicrographPath(md.getValueString(
						MDLabel.MDL_IMAGE, id));
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new TrainingMicrograph(file, ctf, families, getMode());
				loadMicrographData(micrograph);
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No micrographs specified on %s", selfile));
			for (TrainingMicrograph m : micrographs)
				loadAutomaticParticles(m);

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void loadMicrographData(TrainingMicrograph micrograph) {
		try {
			int x, y;
			TrainingParticle particle;
			String filename = micrograph.getOFilename();
			String fname;
			Family family;
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
				if (getMode() == FamilyState.Review
						&& mfd.getStep() != FamilyState.Review) {
					mfd.setState(MicrographFamilyState.Review);
					setChanged(true);
				}
				if (containsBlock(filename, fname)) {
					md2 = new MetaData(fname + "@" + filename);
					for (long id2 : md2.findObjects()) {

						x = md2.getValueInt(MDLabel.MDL_XINT, id2);
						y = md2.getValueInt(MDLabel.MDL_YINT, id2);
						particle = new TrainingParticle(x, y, family, micrograph);
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

	public void persistMicrographs() {
		long id;
		try {
			MetaData md;
			String block = null;
			for (TrainingMicrograph m : micrographs) {
				if (!m.hasData())
					new File(m.getOFilename()).delete();
				else {
					persistMicrographFamilies(m);
					for (MicrographFamilyData mfd : m.getFamiliesData()) {
						if (!mfd.hasManualParticles())
							continue;
						md = new MetaData();
						for (TrainingParticle p : mfd.getManualParticles()) {
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

	public void persistAutomaticParticles(TrainingMicrograph m) {

		if (!m.hasAutomaticParticles())
			new File(getOutputPath(m.getAutoFilename())).delete();
		else
			for (MicrographFamilyData mfd : m.getFamiliesData())
				persistAutomaticParticles(mfd);
	}

	public void loadAutomaticParticles(TrainingMicrograph micrograph) {
		try {
			int x, y;
			double cost;
			boolean deleted = false;
			AutomaticParticle particle;
			MetaData md;
			Family f;
			String file = getOutputPath(micrograph.getAutoFilename());
			if (!new File(file).exists())
				return;
			String[] blocks = MetaData.getBlocksInMetaDataFile(file);
			long[] ids;
			for (String blockname : blocks) {
				md = new MetaData(blockname + "@" + file);
				ids = md.findObjects();
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
				String filename = getOutputPath(mfd.getMicrograph()
						.getOFilename());
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

	public void persistMicrographFamilies(TrainingMicrograph m) {
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

	public void removeFamily(Family family) {
		if (family.particles > 0)
			throw new IllegalArgumentException(
					Constants.getAssociatedDataMsg("family"));
		if (families.size() == 1)
			throw new IllegalArgumentException(
					Constants.getIllegalDeleteMsg("family"));
		families.remove(family);
	}

	public int getNextFreeMicrograph(Family f) {
		int count = 0;
		for (TrainingMicrograph m : micrographs) {
			if (m.getFamilyData(f).getState() == MicrographFamilyState.Available)
				return count;
			count++;
		}
		return -1;
	}



	public void resetFamilyData(MicrographFamilyData mfd) {
		String block;
		MetaData emptymd = new MetaData();
		try {
			// just in case of user reset
			if(this instanceof SupervisedParticlePicker)
				new File(((SupervisedParticlePicker)this).getTrainingAutoFeaturesVectorFile(mfd)).delete();

			if (!mfd.getAutomaticParticles().isEmpty()) {
				// removing automatic particles
				block = String.format("%s@%s", mfd.getFamily().getName(),
						getOutputPath(mfd.getMicrograph().getAutoFilename()));

				emptymd.writeBlock(block);
			}
			if (!mfd.getManualParticles().isEmpty()) {
				// removing manual particles
				block = String.format("%s@%s", mfd.getFamily().getName(), mfd
						.getMicrograph().getOFilename());
				emptymd.writeBlock(block);
			}
			mfd.reset();// Resetting family data

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void saveData() {
		persistFamilies();
		persistMicrographs();

	}

	public int getAutomaticNumber(Family f, double threshold) {
		MicrographFamilyData mfd;
		int count = 0;
		for (TrainingMicrograph m : micrographs) {
			mfd = m.getFamilyData(f);
			count += mfd.getAutomaticParticles(threshold);
		}
		return count;
	}

	
}
