package xmipp.particlepicker.training.model;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;

import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.particlepicker.Family;
import xmipp.particlepicker.Format;
import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePicker;
import xmipp.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.utils.XmippMessage;

public abstract class TrainingPicker extends ParticlePicker {

	protected List<TrainingMicrograph> micrographs;
	private TrainingMicrograph micrograph;

	public static FamilyState previousStep(FamilyState step) {
		if (step == FamilyState.Manual) return null;
		if (step == FamilyState.Supervised) return FamilyState.Manual;
		return null;
	}

	public TrainingPicker(String selfile, String outputdir, String fname, FamilyState mode) {
		super(selfile, outputdir, fname, mode);
		this.micrographs = new ArrayList<TrainingMicrograph>();

	}

	public TrainingPicker(String selfile, String outputdir, FamilyState mode) {
		super(selfile, outputdir, mode);
		this.micrographs = new ArrayList<TrainingMicrograph>();

	}

	public boolean hasEmptyMicrographs(Family f) {
		for (TrainingMicrograph m : micrographs)
			if (m.getFamilyData(f).isEmpty()) return true;
		return false;
	}

	public static FamilyState nextStep(FamilyState step) {
		if (step == FamilyState.Manual) return FamilyState.Supervised;
		if (step == FamilyState.Supervised) return FamilyState.Review;
		return null;
	}

	public List<TrainingMicrograph> getMicrographs() {
		return micrographs;
	}

	public TrainingMicrograph getMicrograph()
	{
		return micrograph;
	}
	

	public void loadMicrographs() {

		micrographs.clear();
		TrainingMicrograph micrograph;
		String ctf = null, file;
		try {
			MetaData md = new MetaData(selfile);
			md.removeDisabled();
			boolean existsctf = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			long[] ids = md.findObjects();
			int fileLabel;

			if (md.containsLabel(MDLabel.MDL_MICROGRAPH))
				fileLabel = MDLabel.MDL_MICROGRAPH;
			else if (md.containsLabel(MDLabel.MDL_IMAGE))
				fileLabel = MDLabel.MDL_IMAGE;
			else
				throw new IllegalArgumentException(String.format("Labels MDL_MICROGRAPH or MDL_IMAGE not found in metadata %s", selfile));

			for (long id : ids) {
				file = md.getValueString(fileLabel, id);
				if (existsctf) ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new TrainingMicrograph(file, ctf, families, getMode());
				loadMicrographData(micrograph);

				micrographs.add(micrograph);
			}
			md.destroy();
			if (micrographs.size() == 0) throw new IllegalArgumentException(String.format("No micrographs specified on %s", selfile));

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void loadMicrographData(TrainingMicrograph micrograph) {
		try {
			String fname;
			Family family;
			MicrographFamilyState state;
			MicrographFamilyData mfd;
			List<MicrographFamilyData> mfdatas = new ArrayList<MicrographFamilyData>();
			if (!new File(getOutputPath(micrograph.getPosFile())).exists()) return;
			MetaData md = new MetaData("families@" + getOutputPath(micrograph.getPosFile()));
			for (long id : md.findObjects()) {

				fname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				state = MicrographFamilyState.valueOf(md.getValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, id));
				family = getFamily(fname);
				mfd = new MicrographFamilyData(micrograph, family, state);
				if (getMode() == FamilyState.Review && mfd.getStep() != FamilyState.Review) {
					mfd.setState(MicrographFamilyState.Review);
					setChanged(true);
				}
				loadManualParticles(mfd, getOutputPath(micrograph.getPosFile()));
				loadAutomaticParticles(mfd, getOutputPath(micrograph.getAutoPosFile()), false);
				mfdatas.add(mfd);
			}
			micrograph.setFamiliesState(mfdatas);
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void loadManualParticles(MicrographFamilyData mfd) {
		loadManualParticles(mfd, getOutputPath(mfd.getMicrograph().getPosFile()));
	}

	public void loadManualParticles(MicrographFamilyData mfd, String file) {
		if (!new File(file).exists()) return;
		Family family = mfd.getFamily();
		if (!containsBlock(file, family.getName())) return;
		int x, y;
		TrainingParticle particle;

		try {
			MetaData md = new MetaData(family.getName() + "@" + file);

			for (long id : md.findObjects()) {

				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				particle = new TrainingParticle(x, y, family, mfd.getMicrograph());
				mfd.addManualParticle(particle);
			}
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void loadAutomaticParticles(MicrographFamilyData mfd) {
		loadAutomaticParticles(mfd, getOutputPath(mfd.getMicrograph().getAutoPosFile()), false);
	}

	public void loadAutomaticParticles(MicrographFamilyData mfd, String file, boolean imported) {
		if (!new File(file).exists()) return;
		Family f = mfd.getFamily();
		if (!containsBlock(file, f.getName())) return;
		int x, y;
		AutomaticParticle particle;
		Double cost;
		boolean deleted;
		try {
			MetaData md = new MetaData(f.getName() + "@" + file);

			for (long id : md.findObjects()) {

				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				cost = md.getValueDouble(MDLabel.MDL_COST, id);
				if (cost == null) throw new IllegalArgumentException("Invalid format for " + file);
				deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? false : true;
				particle = new AutomaticParticle(x, y, f, mfd.getMicrograph(), cost, deleted);
				mfd.addAutomaticParticle(particle, imported);
			}
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void persistMicrographs() {
		try {
			for (TrainingMicrograph m : micrographs) 
				saveData(m);

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void saveData(Micrograph m) {
		TrainingMicrograph tm = (TrainingMicrograph) m;
		long id;
		try {
			MetaData md;
			String block = null;
			String file;
			file = getOutputPath(m.getPosFile());
			if (!m.hasData())
				new File(file).delete();
			else {
				persistMicrographFamilies(tm);
				for (MicrographFamilyData mfd : tm.getFamiliesData()) {
					md = new MetaData();
					for (TrainingParticle p : mfd.getManualParticles()) {
						id = md.addObject();
						md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
						md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
					}
					block = mfd.getFamily().getName() + "@" + file;
					md.writeBlock(block);
					md.destroy();
				}
			}
			persistAutomaticParticles(tm);
			

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void persistAutomaticParticles(TrainingMicrograph m) {

		if (!m.hasAutomaticParticles())
			new File(getOutputPath(m.getAutoPosFile())).delete();
		else
			for (MicrographFamilyData mfd : m.getFamiliesData())
				persistAutomaticParticles(mfd);
	}

	public void persistAutomaticParticles(MicrographFamilyData mfd) {
		try {

			long id;
			if (mfd.hasAutomaticParticles()) {
				String file = getOutputPath(mfd.getMicrograph().getAutoPosFile());
				String section = mfd.getFamily().getName() + "@" + file;
				MetaData md = new MetaData();
				for (AutomaticParticle p : mfd.getAutomaticParticles()) {
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
					md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					md.setValueInt(MDLabel.MDL_ENABLED, (!p.isDeleted()) ? 1 : -1, id);
				}
				md.write(section);
				md.destroy();
			}

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void persistMicrographFamilies(TrainingMicrograph m) {
		long id;
		try {
			String file = getOutputPath(m.getPosFile());
			MetaData md = new MetaData();
			for (MicrographFamilyData mfd : m.getFamiliesData()) {
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, mfd.getFamily().getName(), id);
				md.setValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, mfd.getState().toString(), id);
			}
			md.writeBlock("families@" + file);
			md.destroy();

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public int getNextFreeMicrograph(int index) {
		if (micrographs.size() < index) return -1;
		for (int i = index; i < micrographs.size(); i++) {
			if (micrographs.get(i).getFamilyData(family).getState() == MicrographFamilyState.Available) return i;
		}
		return -1;
	}

	public void resetFamilyData(MicrographFamilyData mfd) {
		String block;
		try {
			MetaData emptymd = new MetaData();
			// just in case of user reset
			if (this instanceof SupervisedParticlePicker)
				new File(((SupervisedParticlePicker) this).getTrainingAutoFeaturesVectorFile(mfd)).delete();

			if (!mfd.getAutomaticParticles().isEmpty()) {
				// removing automatic particles
				block = String.format("%s@%s", mfd.getFamily().getName(), getOutputPath(mfd.getMicrograph().getAutoPosFile()));

				emptymd.writeBlock(block);
			}
			if (!mfd.getManualParticles().isEmpty()) {
				// removing manual particles
				block = String.format("%s@%s", mfd.getFamily().getName(), getOutputPath(mfd.getMicrograph().getPosFile()));
				emptymd.writeBlock(block);
			}
			mfd.reset();// Resetting family data
			emptymd.destroy();

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void saveData() {
		if (isChanged()) {
			super.saveData();
			persistMicrographs();
			// for(Family f: families)
			// {
			// updateFamilyTemplates(f);
			// try {
			// f.getTemplates().write(getOutputPath(f.getName() +
			// "_template.stk"));
			// } catch (Exception e) {
			// getLogger().log(Level.SEVERE, e.getMessage(), e);
			// throw new IllegalArgumentException(e);
			// }
			// }
		}
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

	@Override
	public int getManualParticlesNumber(Family f) {
		int count = 0;
		for (TrainingMicrograph m : micrographs)
			count += m.getFamilyData(f).getManualParticles().size();
		return count;
	}

	public void exportParticles(String file) {

		try {
			MetaData md = new MetaData();
			MicrographFamilyData mfd;
			boolean append = false;
			long id;
			for (TrainingMicrograph m : micrographs) {

				mfd = m.getFamilyData(family);
				if (!mfd.isEmpty()) {
					for (TrainingParticle p : mfd.getParticles()) {
						id = md.addObject();
						md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
						md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
						md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					}
					if (!append)
						md.write("mic_" + m.getName() + "@" + file);
					else
						md.writeBlock("mic_" + m.getName() + "@" + file);
					append = true;
					md.clear();
				}
			}
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	@Override
	public Format detectFormat(String path) {
		Format[] formats = { Format.Xmipp24, Format.Xmipp30, Format.Eman };

		for (TrainingMicrograph m : micrographs) {
			for (Format f : formats) {
				if (Filename.exists(getImportMicrographName(path, m.getFile(), f))) return f;
			}
		}
		return Format.Unknown;
	}

	/** Return the number of particles imported from a file */
	public int importParticlesFromFile(String path, Format f, Micrograph m, float scale, boolean invertx, boolean inverty) {
		MetaData md = new MetaData();
		fillParticlesMdFromFile(path, f, m, md, scale, invertx, inverty);
		int particles = (md != null) ? importParticlesFromMd(m, md) : 0;
		md.destroy();
		return particles;
	}// function importParticlesFromFile

	@Override
	/** Return the number of particles imported */
	public int importParticlesFromFolder(String path, Format f, float scale, boolean invertx, boolean inverty) {
		if (f == Format.Auto) f = detectFormat(path);
		if (f == Format.Unknown) return 0;

		String filename;
		int particles = 0;

		// System.out.println("==========MICROGRAPHS==========");
		// for (TrainingMicrograph m : micrographs)
		// System.out.println("      name: " + m.getFile());
		// System.out.format("  number: %d\n", micrographs.size());
		//
		// System.out.println("==========IMPORTING==========");
		for (TrainingMicrograph m : micrographs) {
			filename = getImportMicrographName(path, m.getFile(), f);
			System.out.println("  filename: " + filename);
			if (Filename.exists(filename)) {
				// System.out.println("    ........EXISTS");
				particles += importParticlesFromFile(filename, f, m, scale, invertx, inverty);
			}
		}
		// System.out.format("==========PARTICLES: %d\n", particles);
		return particles;
	}// function importParticlesFromFolder

	public void importAllParticles(String file) {// Expected a file for all
													// micrographs
		try {
			String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
			List<String> blocks = Arrays.asList(blocksArray);
			String block;
			MetaData md = new MetaData();

			for (TrainingMicrograph m : micrographs) {
				m.reset();
				block = "mic_" + m.getName();
				if (blocks.contains(block)) {
					String blockName = block + "@" + file;
					md.read(blockName);
					importParticlesFromMd(m, md);
				}
			}
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}// function importAllParticles

	/**
	 * Import particles from md, all method to import from files should create
	 * an md and call this function
	 */
	public int importParticlesFromMd(Micrograph m, MetaData md) {
		
		m.reset();
		TrainingMicrograph tm = (TrainingMicrograph) m;
		long[] ids = md.findObjects();
		int x, y;
		double cost;
		boolean hasCost = md.containsLabel(MDLabel.MDL_COST);
		int particles = 0;
		int size = family.getSize();

		for (long id : ids) {
			x = md.getValueInt(MDLabel.MDL_XCOOR, id);
			y = md.getValueInt(MDLabel.MDL_YCOOR, id);
			if (!m.fits(x, y, size))// ignore out of
			// bounds particle
			{
				System.out.println(XmippMessage.getOutOfBoundsMsg("Particle") + String.format(" on x:%s y:%s", x, y));
				continue;
			}
			cost = hasCost ? md.getValueDouble(MDLabel.MDL_COST, id) : 0;
			if (cost == 0 || cost > 1)
				tm.addManualParticle(new TrainingParticle(x, y, family, tm, cost));
			else
				tm.addAutomaticParticle(new AutomaticParticle(x, y, family, tm, cost, false), true);
			++particles;
		}
		return particles;
	}// function importParticlesFromMd
	

	public void removeFamily(Family family) {
		if (getManualParticlesNumber(family) > 0) // perhaps I have to check
													// automatic particles
			throw new IllegalArgumentException(XmippMessage.getAssociatedDataMsg("family"));
		if (families.size() == 1) throw new IllegalArgumentException(XmippMessage.getIllegalDeleteMsg("family"));
		families.remove(family);
		for (TrainingMicrograph m : micrographs)
			m.removeFamilyData(family);
	}

	public static void main(String[] args) {
		try {
			String file = "/home/airen/DNABC/ParticlePicking/Auto/run_001/DefaultFamily_extract_list.xmd";
			MetaData md = new MetaData();
			long[] ids;
			int x, y;
			Double cost;
			String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
			List<String> blocks = Arrays.asList(blocksArray);
			for (String block : blocksArray) {
				String blockName = block + "@" + file;
				md.read(blockName);

				ids = md.findObjects();
				for (long id : ids) {
					x = y = 100;
					cost = 0.0;
					x = md.getValueInt(MDLabel.MDL_XCOOR, id);
					y = md.getValueInt(MDLabel.MDL_YCOOR, id);
					cost = md.getValueDouble(MDLabel.MDL_COST, id);
				}
			}
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	// public void updateFamilyTemplates(Family f) {
	// ImageGeneric igp;
	// List<TrainingParticle> particles;
	// MicrographFamilyData mfd;
	// for(TrainingMicrograph m: micrographs)
	// {
	// mfd = m.getFamilyData(f);
	// for (int i = 0; i < mfd.getManualParticles().size(); i++) {
	// particles = mfd.getManualParticles();
	// igp = particles.get(i).getImageGeneric();
	// if (i < f.getTemplatesNumber())
	// f.setTemplate((int) (ImageGeneric.FIRST_IMAGE + i), igp);
	// else
	// try {
	// f.getTemplates().alignImages(igp);
	// } catch (Exception e) {
	// throw new IllegalArgumentException(e.getMessage());
	// }
	// }
	// }
	//
	// }
	
	public String getImportMicrographName(String path, String filename, Format f) {
		String base = Filename.removeExtension(Filename.getBaseName(filename));
		switch (f) {
		case Xmipp24:
			return Filename.join(path, base, base + ".raw.Common.pos");
		case Xmipp30:
			return Filename.join(path, base + ".pos");
		case Eman:
			return Filename.join(path, base + ".box");

		default:
			return null;
		}
	}
	
	@Override
	public void setMicrograph(Micrograph m) {
		this.micrograph = (TrainingMicrograph)m;
		
	}

}
