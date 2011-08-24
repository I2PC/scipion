package model;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import xmipp.MDLabel;
import xmipp.MetaData;
import xmipp.Program;






public class ParticlePicker {

	private List<Family> families;
	private List<Micrograph> micrographs;
	private static ParticlePicker ppicker;

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
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
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
				md.setValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, f.getStep()
						.toString(), id);
			}
			md.write(filename);
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
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
						MDLabel.MDL_ASSOCIATED_IMAGE1, id));
				family = new Family(gname, new Color(rgb), size, step);
				families.add(family);
			}
			if (families.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No families specified on %s", filename));
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
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
		String ctf, filename;
		try {
			MetaData md = new MetaData(xmd);
			long[] ids = md.findObjects();
			for (long id : ids) {

				filename = ExecutionEnvironment.getMicrographPath(md
						.getValueString(MDLabel.MDL_IMAGE, id));
				ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new Micrograph(filename, ctf);
				loadMicrographData(micrograph);
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No micrographs specified on %s", xmd));

		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
					e);
		}

	}

	public void loadMicrographData(Micrograph micrograph) {
		try {
			int x, y;
			Particle particle;
			String filename = micrograph.getOutputFName();
			String fname, step;
			Family family;
			MicrographFamilyData mfd;
			List<MicrographFamilyData> mfdatas = new ArrayList<MicrographFamilyData>();
			if (!new File(filename).exists())
				return;
			MetaData md = new MetaData("families@" + filename);
			MetaData md2;
			for (long id : md.findObjects()) {
				fname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				step = md.getValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, id);
				family = getFamily(fname);
				mfd = new MicrographFamilyData(family, Step.valueOf(step));
				md2 = new MetaData(fname + "@" + filename);
				for (long id2 : md2.findObjects()) {

					x = md2.getValueInt(MDLabel.MDL_XINT, id2);
					y = md2.getValueInt(MDLabel.MDL_YINT, id2);
					particle = new Particle(x, y, family, micrograph);
					mfd.addParticle(particle);
				}
				mfdatas.add(mfd);
			}
			micrograph.setFamiliesData(mfdatas);
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
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
					md.setValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, mfd
							.getStep().toString(), id);
				}
				md.write("families@" + m.getOutputFName());
				for (MicrographFamilyData mfd : m.getFamiliesData()) {
					if (mfd.isEmpty())
						continue;
					md = new MetaData();
					for (Particle p : mfd.getParticles()) {
						id = md.addObject();
						md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
						md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
					}
					block = mfd.getFamily().getName() + "@"
							+ m.getOutputFName();
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
			Particle particle;
			MetaData md;
			Family f;
			String filename = micrograph.getAutoOutputFName();
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
					particle = new Particle(x, y, f, micrograph, true, cost);
					micrograph.addParticle(particle);
				}
			}
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
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
								ExecutionEnvironment.getOutputPath(f.getName()),// --model
								ExecutionEnvironment.getOutputPath(micrograph
										.getName()), // --outputRoot
								f.getName() + "@" + micrograph.getOutputFName());// train
																					// parameter
				if (ExecutionEnvironment.isFastMode())
					args += " --fast";
				if (ExecutionEnvironment.isIncore())
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
						ExecutionEnvironment.getOutputPath(f.getName()),// --model
						ExecutionEnvironment
								.getOutputPath(micrograph.getName()),// --outputRoot
						ExecutionEnvironment.getThreads()// --thr
				);

		if (ExecutionEnvironment.isFastMode())
			args += " --fast";
		if (ExecutionEnvironment.isIncore())
			args += " --in_core";
		executeProgram("xmipp_micrograph_automatic_picking", args);
	}

}
