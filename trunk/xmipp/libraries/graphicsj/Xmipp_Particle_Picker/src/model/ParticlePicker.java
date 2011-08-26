package model;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
	private static int mintraining = 20;
	private static String trainingfn = "training";
	private static String trainingmaskfn = "mask";
	private static double mincost = 0;
	
	private List<Family> families;
	private List<Micrograph> micrographs;
	private static ParticlePicker ppicker;
	
	public static Logger getLogger()
	{
		try {
			if(logger == null)
			{
				FileHandler fh = new FileHandler("PPicker.log", true);
				fh.setFormatter(new SimpleFormatter());
				logger = Logger.getLogger("PPickerLogger");
				logger.addHandler(fh);
			}
			return logger;
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	public static String getTrainingFilename()
	{
		return trainingfn;
	}
	
	public static String getTrainingMaskFilename()
	{
		return trainingmaskfn;
	}
	
	public static int getMinForTraining()
	{
		return mintraining;
	}
	
	public static boolean getIsAuto()
	{
		return auto;
	}
	
	public static int getMinParticles()
	{
		return minparticles;
	}
	
	public static String getOutputDir()
	{
		return outputdir;
	}
	
	
	public static String getMicrographPath(String rpath)
	{
		return rundir + File.separator + rpath;
	}
	
	public static void setMicrographsSelFile(String file)
	{
		mgselfile = file;
	}
	
	public static String getMicrographsSelFile()
	{
		return mgselfile;
	}
	
	public static void setOutputDir(String dir)
	{
		if(!new File(dir).exists())
			throw new IllegalArgumentException("Output dir " + dir + " does not exist");
		outputdir = dir;
	}
	
	public static String getOutputPath(String file)
	{
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
	
	public static boolean getFastMode()
	{
		return fastMode;
	}

	public static void setIsAuto(boolean automatic) {
		ParticlePicker.auto = automatic;
		
	}
	
	public static int getThreads()
	{
		return threads;
	}


	public boolean hasEmptyMicrographs(Family f)
	{
		for(Micrograph m: micrographs)
			if(m.getFamilyData(f).isEmpty())
				return true;
		return false;
	}

	private ParticlePicker() {
		this.families = new ArrayList<Family>();
		this.micrographs = new ArrayList<Micrograph>();
		loadFamilyData();
		loadMicrographsData();
	}

	public static Step nextStep(Step step) {
		if (step == Step.Manual)
			return Step.Supervised;
		if (step == Step.Supervised)
			return Step.Manual;
		return null;
	}

	private void executeProgram(String name, String args) {
		System.out.println(args);
		try {
			Program.runByName(name, args);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(),
					e);
			throw new IllegalArgumentException(e.getMessage());
		}
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
			getLogger().log(Level.SEVERE, e.getMessage(),
					e);
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
		Step step;
		String gname;
		try {
			MetaData md = new MetaData(filename);
			long[] ids = md.findObjects();
			for (long id : ids) {
				gname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				rgb = md.getValueInt(MDLabel.MDL_PICKING_COLOR, id);
				size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
				step = Step.valueOf(md.getValueString(
						MDLabel.MDL_PICKING_FAMILY_STATE, id));
				family = new Family(gname, new Color(rgb), size, step);
				families.add(family);
			}
			if (families.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No families specified on %s", filename));
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(),
					e);
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

				filename = getMicrographPath(md
						.getValueString(MDLabel.MDL_IMAGE, id));
				System.out.println(filename);
				if(existsctf)
					ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new Micrograph(filename, ctf);
				loadMicrographData(micrograph);
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No micrographs specified on %s", xmd));

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(),
					e);
		}

	}

	public void loadMicrographData(Micrograph micrograph) {
		try {
			int x, y;
			Particle particle;
			String filename = micrograph.getOFilename();
			String fname;
			Family family;
			Step step;
			State state;
			MicrographFamilyData mfd;
			List<MicrographFamilyData> mfdatas = new ArrayList<MicrographFamilyData>();
			if (!new File(filename).exists())
				return;
			MetaData md = new MetaData("families@" + filename);
			MetaData md2;
			for (long id : md.findObjects()) {
				fname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				step = Step.valueOf(md.getValueString(MDLabel.MDL_PICKING_FAMILY_STATE, id));
				state = State.valueOf(md.getValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, id));
				family = getFamily(fname);
				mfd = new MicrographFamilyData(micrograph, family, step, state);
				md2 = new MetaData(fname + "@" + filename);
				for (long id2 : md2.findObjects()) {

					x = md2.getValueInt(MDLabel.MDL_XINT, id2);
					y = md2.getValueInt(MDLabel.MDL_YINT, id2);
					particle = new Particle(x, y, family, micrograph);
					mfd.addManualParticle(particle);
				}
				mfdatas.add(mfd);
			}
			micrograph.setFamiliesData(mfdatas);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(),
					e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void persistMicrographsData() {
		long id;
		try {
			MetaData md;
			String block = null;
			for (Micrograph m : micrographs) {
				if (m.isEmpty())
					continue;
				md = new MetaData();
				for (MicrographFamilyData mfd : m.getFamiliesData()) {
					if (mfd.isEmpty())
						continue;
					id = md.addObject();
					md.setValueString(MDLabel.MDL_PICKING_FAMILY, mfd
							.getFamily().getName(), id);
					md.setValueString(MDLabel.MDL_PICKING_FAMILY_STATE, mfd
							.getStep().toString(), id);
					md.setValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, mfd.getState().toString(), id);
				}
				md.write("families@" + m.getOFilename());
				for (MicrographFamilyData mfd : m.getFamiliesData()) {
					if (mfd.isEmpty())
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

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
			System.out.println(filename);
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
					deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1)? false:true;					
					particle = new AutomaticParticle(x, y, f, micrograph, cost, deleted);
					micrograph.addAutomaticParticle(particle);
				}
			}
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(),
					e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}
	
	public void persistAutomaticParticles(MicrographFamilyData mfd)
	{
		try {
			long id;
			String filename = mfd.getMicrograph().getAutoOFilename();

			MetaData md = new MetaData();
			for (AutomaticParticle p : mfd.getAutomaticParticles()) {
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
					md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					md.setValueInt(MDLabel.MDL_ENABLED, (!p.isDeleted())? 1: -1, id);
			}
			md.write(mfd.getFamily().getName() + "@" + filename);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(),
					e);
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
		if(micrograph.getFamilyData(family).getManualParticles().size() > 0)
			args += family.getName() + "@" + micrograph.getOFilename();
		if (isFastMode())
			args += " --fast";
		if (isIncore())
			args += " --in_core";
		executeProgram("xmipp_micrograph_automatic_picking", args);
		
	}

}
