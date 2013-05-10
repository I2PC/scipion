package xmipp.viewer.particlepicker.training.model;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.jni.Particle;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Family;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;

public abstract class TrainingPicker extends ParticlePicker
{

	protected List<TrainingMicrograph> micrographs;
	private TrainingMicrograph micrograph;

	public static FamilyState previousStep(FamilyState step)
	{
		if (step == FamilyState.Manual)
			return null;
		if (step == FamilyState.Supervised)
			return FamilyState.Manual;
		return null;
	}

	public TrainingPicker(String selfile, String outputdir, String fname, FamilyState mode)
	{
		super(selfile, outputdir, fname, mode);

	}

	public TrainingPicker(String selfile, String outputdir, FamilyState mode)
	{
		super(selfile, outputdir, mode);

	}

	public boolean hasEmptyMicrographs(Family f)
	{
		for (TrainingMicrograph m : micrographs)
			if (m.getFamilyData(f).isEmpty())
				return true;
		return false;
	}

	public static FamilyState nextStep(FamilyState step)
	{
		if (step == FamilyState.Manual)
			return FamilyState.Supervised;
		if (step == FamilyState.Supervised)
			return FamilyState.Review;
		return null;
	}

	public List<TrainingMicrograph> getMicrographs()
	{
		return micrographs;
	}

	public TrainingMicrograph getMicrograph()
	{
		return micrograph;
	}

	public void loadMicrographs()
	{
		try
		{
			loadEmptyMicrographs();
			for (TrainingMicrograph m : micrographs)
				loadMicrographData(m);

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void loadMicrographData(TrainingMicrograph micrograph)
	{
		try
		{
			String fname;
			Family family;
			MicrographFamilyState state;
			MicrographFamilyData mfd;
			Integer autopickpercent;
			List<MicrographFamilyData> mfdatas = new ArrayList<MicrographFamilyData>();
			if (!new File(getOutputPath(micrograph.getPosFile())).exists())
				return;

			MetaData md = new MetaData("families@" + getOutputPath(micrograph.getPosFile()));
			boolean hasautopercent = md.containsLabel(MDLabel.MDL_PICKING_AUTOPICKPERCENT);
			for (long id : md.findObjects())
			{

				fname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				state = MicrographFamilyState.valueOf(md.getValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, id));
				family = getFamily(fname);
				if (family == null)
					throw new IllegalArgumentException(XmippMessage.getIllegalValueMsg("family", fname));
				if (hasautopercent)
					autopickpercent = md.getValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, id);
				else
					autopickpercent = 50;// compatibility with previous projects
				mfd = new MicrographFamilyData(micrograph, family, state, autopickpercent);

				if (getMode() == FamilyState.Review && mfd.getStep() != FamilyState.Review)
				{
					mfd.setState(MicrographFamilyState.Review);
					setChanged(true);
				}
				loadManualParticles(mfd, getOutputPath(micrograph.getPosFile()));
				loadAutomaticParticles(mfd, getOutputPath(micrograph.getAutoPosFile()), false);
				mfdatas.add(mfd);
			}
			micrograph.setFamiliesState(mfdatas);
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void loadManualParticles(MicrographFamilyData mfd)
	{
		loadManualParticles(mfd, getOutputPath(mfd.getMicrograph().getPosFile()));
	}

	public void loadManualParticles(MicrographFamilyData mfd, String file)
	{
		if (!new File(file).exists())
			return;
		Family family = mfd.getFamily();
		if (!containsBlock(file, family.getName()))
			return;
		int x, y;
		TrainingParticle particle;

		try
		{
			MetaData md = new MetaData(family.getName() + "@" + file);

			for (long id : md.findObjects())
			{

				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				particle = new TrainingParticle(x, y, family, mfd.getMicrograph());
				mfd.addManualParticle(particle, this, false, false);
			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void loadAutomaticParticles(MicrographFamilyData mfd)
	{
		loadAutomaticParticles(mfd, getOutputPath(mfd.getMicrograph().getAutoPosFile()), false);
	}

	public void loadAutomaticParticles(MicrographFamilyData mfd, String file, boolean imported)
	{
		if (!new File(file).exists())
			return;
		Family f = mfd.getFamily();
		if (!containsBlock(file, f.getName()))
			return;
		int x, y;
		AutomaticParticle particle;
		Double cost;
		boolean deleted;
		try
		{
			MetaData md = new MetaData(f.getName() + "@" + file);

			for (long id : md.findObjects())
			{

				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				cost = md.getValueDouble(MDLabel.MDL_COST, id);
				if (cost == null)
					throw new IllegalArgumentException("Invalid format for " + file);
				deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? false : true;
				particle = new AutomaticParticle(x, y, f, mfd.getMicrograph(), cost, deleted);
				mfd.addAutomaticParticle(particle, imported);
			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void saveMicrographs()
	{
		try
		{
			for (TrainingMicrograph m : micrographs)
				saveData(m);

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void saveData(Micrograph m)
	{
		TrainingMicrograph tm = (TrainingMicrograph) m;
		long id;
		try
		{
			MetaData md;
			String block = null;
			String file;
			file = getOutputPath(m.getPosFile());
			if (!m.hasData())
				new File(file).delete();
			else
			{
				persistMicrographState(tm);
				for (MicrographFamilyData mfd : tm.getFamiliesData())
				{
					md = new MetaData();
					for (TrainingParticle p : mfd.getManualParticles())
					{
						id = md.addObject();
						md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
						md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
					}
					block = mfd.getFamily().getName() + "@" + file;
					md.writeBlock(block);
					md.destroy();
				}
			}
			saveAutomaticParticles(tm);
			// saveTemplates();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void saveAutomaticParticles(TrainingMicrograph m)
	{

		if (!m.hasAutomaticParticles())
			new File(getOutputPath(m.getAutoPosFile())).delete();
		else
			for (MicrographFamilyData mfd : m.getFamiliesData())
				saveAutomaticParticles(mfd);
	}

	public void saveAutomaticParticles(MicrographFamilyData mfd)
	{
		try
		{

			long id;
			if (mfd.hasAutomaticParticles())
			{
				String file = getOutputPath(mfd.getMicrograph().getAutoPosFile());
				String section = mfd.getFamily().getName() + "@" + file;
				MetaData md = new MetaData();
				for (AutomaticParticle p : mfd.getAutomaticParticles())
				{
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
					md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					md.setValueInt(MDLabel.MDL_ENABLED, (!p.isDeleted()) ? 1 : -1, id);
				}
				md.write(section);
				md.destroy();
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void persistMicrographState(TrainingMicrograph m)
	{
		long id;
		try
		{
			String file = getOutputPath(m.getPosFile());
			MetaData md = new MetaData();
			for (MicrographFamilyData mfd : m.getFamiliesData())
			{
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, mfd.getFamily().getName(), id);
				md.setValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, mfd.getState().toString(), id);
				md.setValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, mfd.getAutopickpercent(), id);
			}
			md.writeBlock("families@" + file);
			md.destroy();

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public int getNextFreeMicrograph(int index)
	{
		if (micrographs.size() < index)
			return -1;
		for (int i = index; i < micrographs.size(); i++)
		{
			if (micrographs.get(i).getFamilyData(family).getState() == MicrographFamilyState.Available)
				return i;
		}
		return -1;
	}

	public void resetFamilyData(MicrographFamilyData mfd)
	{
		String block;
		try
		{
			MetaData emptymd = new MetaData();
			// just in case of user reset
			if (this instanceof SupervisedParticlePicker)
				new File(((SupervisedParticlePicker) this).getTrainingAutoFeaturesVectorFile(mfd)).delete();

			if (!mfd.getAutomaticParticles().isEmpty())
			{
				// removing automatic particles
				block = String.format("%s@%s", mfd.getFamily().getName(), getOutputPath(mfd.getMicrograph().getAutoPosFile()));

				emptymd.writeBlock(block);
			}
			if (!mfd.getManualParticles().isEmpty())
			{
				// removing manual particles
				block = String.format("%s@%s", mfd.getFamily().getName(), getOutputPath(mfd.getMicrograph().getPosFile()));
				emptymd.writeBlock(block);
			}
			mfd.reset();// Resetting family data
			emptymd.destroy();

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void saveData()
	{
		if (isChanged())
		{
			super.saveData();
			saveMicrographs();
		}
		if (getMode() == FamilyState.Manual)// only changed in manual mode
			saveTemplates();
	}

	public int getAutomaticNumber(Family f, double threshold)
	{
		MicrographFamilyData mfd;
		int count = 0;
		for (TrainingMicrograph m : micrographs)
		{
			mfd = m.getFamilyData(f);
			count += mfd.getAutomaticParticles(threshold);
		}
		return count;
	}

	@Override
	public int getManualParticlesNumber(Family f)
	{
		int count = 0;
		for (TrainingMicrograph m : micrographs)
			count += m.getFamilyData(f).getManualParticles().size();
		return count;
	}

	public void exportParticles(String file)
	{

		try
		{
			MetaData md = new MetaData();
			MicrographFamilyData mfd;
			boolean append = false;
			long id;
			for (TrainingMicrograph m : micrographs)
			{

				mfd = m.getFamilyData(family);
				if (!mfd.isEmpty())
				{
					for (TrainingParticle p : mfd.getParticles())
					{
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
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	@Override
	public Format detectFormat(String path)
	{
		Format[] formats = { Format.Xmipp24, Format.Xmipp30, Format.Eman };

		for (TrainingMicrograph m : micrographs)
		{
			for (Format f : formats)
			{
				if (Filename.exists(getImportMicrographName(path, m.getFile(), f)))
					return f;
			}
		}
		return Format.Unknown;
	}

	/** Return the number of particles imported from a file */
	public String importParticlesFromFile(String path, Format f, Micrograph m, float scale, boolean invertx, boolean inverty)
	{
		MetaData md = new MetaData();
		fillParticlesMdFromFile(path, f, m, md, scale, invertx, inverty);
		String result = (md != null) ? importParticlesFromMd(m, md) : "";
		md.destroy();
		return result;
	}// function importParticlesFromFile

	@Override
	/** Return the number of particles imported */
	public String importParticlesFromFolder(String path, Format f, float scale, boolean invertx, boolean inverty)
	{
		if (f == Format.Auto)
			f = detectFormat(path);
		if (f == Format.Unknown)
			return "Unknown format";

		String filename;
		String result = "";
		for (TrainingMicrograph m : micrographs)
		{
			filename = getImportMicrographName(path, m.getFile(), f);
			System.out.println("  filename: " + filename);
			if (Filename.exists(filename))
				result += importParticlesFromFile(filename, f, m, scale, invertx, inverty);
		}
		return result;
	}// function importParticlesFromFolder

	public void importAllParticles(String file)
	{// Expected a file for all
		// micrographs
		try
		{
			String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
			List<String> blocks = Arrays.asList(blocksArray);
			String block;
			MetaData md = new MetaData();

			for (TrainingMicrograph m : micrographs)
			{
				m.reset();
				block = "mic_" + m.getName();
				if (blocks.contains(block))
				{
					String blockName = block + "@" + file;
					md.read(blockName);
					importParticlesFromMd(m, md);
				}
			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}// function importAllParticles

	public String importAllParticles(String file, float scale, boolean invertx, boolean inverty)
	{// Expected a file for all
		// micrographs
		try
		{
			String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
			List<String> blocks = Arrays.asList(blocksArray);
			String block;
			MetaData md = new MetaData();
			int width, height;

			for (TrainingMicrograph m : micrographs)
			{
				m.reset();
				block = "mic_" + m.getName();
				if (blocks.contains(block))
				{
					String blockName = block + "@" + file;
					md.read(blockName);
					width = (int) (m.width / scale);// original width
					height = (int) (m.height / scale);// original height
					if (invertx)
						md.operate(String.format("xcoor=%d-xcoor", width));
					if (inverty)
						md.operate(String.format("ycoor=%d-ycoor", height));
					if (scale != 1.f)
						md.operate(String.format("xcoor=xcoor*%f,ycoor=ycoor*%f", scale, scale));
					importParticlesFromMd(m, md);
				}
			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
		return null;

	}// function importAllParticles

	/**
	 * Import particles from md, all method to import from files should create
	 * an md and call this function
	 */
	public String importParticlesFromMd(Micrograph m, MetaData md)
	{

		m.reset();
		TrainingMicrograph tm = (TrainingMicrograph) m;
		long[] ids = md.findObjects();
		int x, y;
		double cost;
		boolean hasCost = md.containsLabel(MDLabel.MDL_COST);
		String result = "";
		int size = family.getSize();

		for (long id : ids)
		{
			x = md.getValueInt(MDLabel.MDL_XCOOR, id);
			y = md.getValueInt(MDLabel.MDL_YCOOR, id);
			if (!m.fits(x, y, size))// ignore out of
			// bounds particle
			{
				result += XmippMessage.getOutOfBoundsMsg("Particle") + String.format(" on x:%s y:%s", x, y);
				continue;
			}
			cost = hasCost ? md.getValueDouble(MDLabel.MDL_COST, id) : 0;
			if (cost == 0 || cost > 1)
				tm.addManualParticle(new TrainingParticle(x, y, family, tm, cost), this, false, true);
			else
				tm.addAutomaticParticle(new AutomaticParticle(x, y, family, tm, cost, false), true);
		}
		return result;
	}// function importParticlesFromMd

	public void removeFamily(Family family)
	{
		if (getManualParticlesNumber(family) > 0) // perhaps I have to check
													// automatic particles
			throw new IllegalArgumentException(XmippMessage.getAssociatedDataMsg("family"));
		if (families.size() == 1)
			throw new IllegalArgumentException(XmippMessage.getIllegalDeleteMsg("family"));
		families.remove(family);
		for (TrainingMicrograph m : micrographs)
			m.removeFamilyData(family);
	}

	public static void main(String[] args)
	{
		try
		{
			String file = "/home/airen/DNABC/ParticlePicking/Auto/run_001/DefaultFamily_extract_list.xmd";
			MetaData md = new MetaData();
			long[] ids;
			int x, y;
			Double cost;
			String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
			List<String> blocks = Arrays.asList(blocksArray);
			for (String block : blocksArray)
			{
				String blockName = block + "@" + file;
				md.read(blockName);

				ids = md.findObjects();
				for (long id : ids)
				{
					x = y = 100;
					cost = 0.0;
					x = md.getValueInt(MDLabel.MDL_XCOOR, id);
					y = md.getValueInt(MDLabel.MDL_YCOOR, id);
					cost = md.getValueDouble(MDLabel.MDL_COST, id);
				}
			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	public String getImportMicrographName(String path, String filename, Format f)
	{
		String base = Filename.removeExtension(Filename.getBaseName(filename));
		switch (f)
		{
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
	public void setMicrograph(Micrograph m)
	{
		this.micrograph = (TrainingMicrograph) m;

	}

	public void loadEmptyMicrographs()
	{
		if (micrographs == null)
			micrographs = new ArrayList<TrainingMicrograph>();
		else
			micrographs.clear();
		TrainingMicrograph micrograph;
		String psd = null, ctf = null, filename;
		try
		{
			MetaData md = new MetaData(getMicrographsSelFile());
			md.removeDisabled();
			int fileLabel;

			if (md.containsLabel(MDLabel.MDL_MICROGRAPH))
				fileLabel = MDLabel.MDL_MICROGRAPH;
			else if (md.containsLabel(MDLabel.MDL_IMAGE))
				fileLabel = MDLabel.MDL_IMAGE;
			else
				throw new IllegalArgumentException(String.format("Labels MDL_MICROGRAPH or MDL_IMAGE not found in metadata %s", selfile));
			boolean existspsd = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			boolean existsctf = md.containsLabel(MDLabel.MDL_CTF_MODEL);
			long[] ids = md.findObjects();
			for (long id : ids)
			{

				filename = md.getValueString(fileLabel, id);
				if (existspsd)
					psd = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_CTF_MODEL, id);
				micrograph = new TrainingMicrograph(filename, psd, ctf, families, getMode());
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", getMicrographsSelFile()));
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	public MicrographFamilyData getFamilyData()
	{
		return micrograph.getFamilyData(family);
	}

	public boolean isReviewFile(String file)
	{
		try
		{
			MetaData md = new MetaData(file);

			if (!md.containsLabel(MDLabel.MDL_XCOOR))
				return false;
			if (!md.containsLabel(MDLabel.MDL_YCOOR))
				return false;
			if (!md.containsLabel(MDLabel.MDL_COST))
				return false;
			String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
			List<String> blocks = Arrays.asList(blocksArray);

			String name;
			for (String block : blocks)
			{
				name = block.replace("mic_", "");
				if (getMicrograph(name) == null)
					return false;

			}

			md.destroy();
			return true;
		}
		catch (Exception e)
		{
			return false;
		}
	}

	public void updateTemplates()
	{
		updateTemplates(family);

	}

	public void updateTemplates(Family f)
	{
		if (f.getStep() != FamilyState.Manual)
			return;

		if (!hasManualParticles(f))
			return;

		f.initTemplates();
		ImageGeneric igp;
		List<TrainingParticle> particles;
		MicrographFamilyData mfd;
		TrainingParticle particle;
		try
		{
			for (TrainingMicrograph m : micrographs)
			{
				mfd = m.getFamilyData(f);
				for (int i = 0; i < mfd.getManualParticles().size(); i++)
				{
					particles = mfd.getManualParticles();
					particle = particles.get(i);
					igp = particle.getImageGeneric();
					if (f.getTemplateIndex() < f.getTemplatesNumber())
						f.setTemplate(igp);
					else
					{
						double[] align = family.getTemplates().alignImage(igp);
						particle.setLastalign(align);
					}
				}
			}
			f.saveTemplates();
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void saveTemplates()
	{
		try
		{
			for (Family f : families)
				f.saveTemplates();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void resetParticleImages()
	{
		MicrographFamilyData mfd;
		for (TrainingMicrograph m : micrographs)
		{
			mfd = m.getFamilyData(family);
			for (TrainingParticle p : mfd.getManualParticles())
				p.resetImagePlus();

		}
	}

	public void loadConfig()
	{
		String file = configfile;
		if (!new File(file).exists())
		{
			if (family == null)
				family = families.get(0);
			setMicrograph(getMicrographs().get(0));
			return;

		}

		String mname, fname;
		try
		{
			MetaData md = new MetaData(file);
			boolean hasautopercent = md.containsLabel(MDLabel.MDL_PICKING_AUTOPICKPERCENT);
			for (long id : md.findObjects())
			{
				if (family == null)
				{
					fname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
					family = getFamily(fname);
				}
				mname = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
				setMicrograph(getMicrograph(mname));
				if (hasautopercent)
					autopickpercent = md.getValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, id);

			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void saveConfig()
	{
		try
		{
			MetaData md;
			String file = configfile;
			md = new MetaData();
			long id = md.addObject();
			md.setValueString(MDLabel.MDL_PICKING_FAMILY, family.getName(), id);
			md.setValueString(MDLabel.MDL_MICROGRAPH, getMicrograph().getName(), id);
			md.setValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, getAutopickpercent(), id);
			md.write(file);
			md.destroy();

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public boolean hasManualParticles()
	{
		return hasManualParticles(family);
	}

	public boolean hasManualParticles(Family f)
	{
		MicrographFamilyData mfd = null;
		for (TrainingMicrograph m : micrographs)
		{
			mfd = m.getFamilyData(f);
			if (mfd.hasManualParticles())
				return true;
		}
		return false;
	}

	@Override
	public boolean isValidSize(int size)
	{
		for (TrainingParticle p : getFamilyData().getParticles())
			if (!getMicrograph().fits(p.getX(), p.getY(), size))
				return false;
		return true;
	}

}
