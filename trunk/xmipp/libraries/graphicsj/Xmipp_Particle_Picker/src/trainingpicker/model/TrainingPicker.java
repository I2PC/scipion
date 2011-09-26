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

public abstract class TrainingPicker extends ParticlePicker
{

	private String selfile = "micrographs.sel";

	protected List<TrainingMicrograph> micrographs;

	public static FamilyState previousStep(FamilyState step)
	{
		if (step == FamilyState.Manual)
			return null;
		if (step == FamilyState.Supervised)
			return FamilyState.Manual;
		return null;
	}

	public TrainingPicker(String selfile, String outputdir, FamilyState mode)
	{
		super(outputdir, mode);
		this.selfile = selfile;
		this.micrographs = new ArrayList<TrainingMicrograph>();

	}

	public String getMicrographsSelFile()
	{
		return selfile;
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

	public TrainingMicrograph getMicrograph(String name)
	{
		for (TrainingMicrograph m : getMicrographs())
			if (m.getName().equalsIgnoreCase(name))
				return m;
		return null;
	}

	public void loadMicrographs()
	{

		micrographs.clear();
		TrainingMicrograph micrograph;
		String ctf = null, file;
		try
		{
			MetaData md = new MetaData(selfile);
			boolean existsctf = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			long[] ids = md.findObjects();
			for (long id : ids)
			{

				file = md.getValueString(MDLabel.MDL_IMAGE, id);
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new TrainingMicrograph(file, ctf, families, getMode());
				loadMicrographData(micrograph);
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", selfile));
			for (TrainingMicrograph m : micrographs)
				loadAutomaticParticles(m);

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
			int x, y;
			TrainingParticle particle;
			String file = getOutputPath(micrograph.getOFilename());
			String fname;
			Family family;
			MicrographFamilyState state;
			MicrographFamilyData mfd;
			List<MicrographFamilyData> mfdatas = new ArrayList<MicrographFamilyData>();
			if (!new File(file).exists())
				return;
			MetaData md = new MetaData("families@" + file);
			MetaData md2;
			for (long id : md.findObjects())
			{
				fname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				state = MicrographFamilyState.valueOf(md.getValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, id));
				family = getFamily(fname);
				mfd = new MicrographFamilyData(micrograph, family, state);
				if (getMode() == FamilyState.Review && mfd.getStep() != FamilyState.Review)
				{
					mfd.setState(MicrographFamilyState.Review);
					setChanged(true);
				}
				if (containsBlock(file, fname))
				{
					md2 = new MetaData(fname + "@" + file);
					for (long id2 : md2.findObjects())
					{

						x = md2.getValueInt(MDLabel.MDL_XINT, id2);
						y = md2.getValueInt(MDLabel.MDL_YINT, id2);
						particle = new TrainingParticle(x, y, family, micrograph);
						mfd.addManualParticle(particle);
					}
				}

				mfdatas.add(mfd);
			}
			micrograph.setFamiliesState(mfdatas);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void persistMicrographs()
	{
		long id;
		try
		{
			MetaData md;
			String block = null;
			for (TrainingMicrograph m : micrographs)
			{
				if (!m.hasData())
					new File(getOutputPath(m.getOFilename())).delete();
				else
				{
					persistMicrographFamilies(m);
					for (MicrographFamilyData mfd : m.getFamiliesData())
					{
						if (!mfd.hasManualParticles())
							continue;
						md = new MetaData();
						for (TrainingParticle p : mfd.getManualParticles())
						{
							id = md.addObject();
							md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
							md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
						}
						block = mfd.getFamily().getName() + "@" + m.getOFilename();
						md.writeBlock(block);
					}
				}
				persistAutomaticParticles(m);
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void persistAutomaticParticles(TrainingMicrograph m)
	{

		if (!m.hasAutomaticParticles())
			new File(getOutputPath(m.getAutoFilename())).delete();
		else
			for (MicrographFamilyData mfd : m.getFamiliesData())
				persistAutomaticParticles(mfd);
	}

	public void loadAutomaticParticles(TrainingMicrograph micrograph)
	{
		try
		{
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
			for (String blockname : blocks)
			{
				md = new MetaData(blockname + "@" + file);
				ids = md.findObjects();
				for (long id : ids)
				{

					x = md.getValueInt(MDLabel.MDL_XINT, id);
					y = md.getValueInt(MDLabel.MDL_YINT, id);
					cost = md.getValueDouble(MDLabel.MDL_COST, id);
					f = getFamily(blockname);
					if (f == null)
						throw new IllegalArgumentException("Unknown family " + blockname);
					deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? false : true;
					particle = new AutomaticParticle(x, y, f, micrograph, cost, deleted);
					micrograph.addAutomaticParticle(particle);
				}
			}
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void persistAutomaticParticles(MicrographFamilyData mfd)
	{
		try
		{
			long id;
			if (mfd.getMicrograph().hasAutomaticParticles())
			{
				String filename = getOutputPath(mfd.getMicrograph().getOFilename());
				MetaData md = new MetaData();
				for (AutomaticParticle p : mfd.getAutomaticParticles())
				{
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
					md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					md.setValueInt(MDLabel.MDL_ENABLED, (!p.isDeleted()) ? 1 : -1, id);
				}
				md.write(mfd.getFamily().getName() + "@" + filename);
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void persistMicrographFamilies(TrainingMicrograph m)
	{
		long id;
		try
		{
			MetaData md = new MetaData();
			for (MicrographFamilyData mfd : m.getFamiliesData())
			{
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, mfd.getFamily().getName(), id);
				md.setValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, mfd.getState().toString(), id);
			}
			md.writeBlock("families@" + m.getOFilename());

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public int getNextFreeMicrograph(Family f)
	{
		int count = 0;
		for (TrainingMicrograph m : micrographs)
		{
			if (m.getFamilyData(f).getState() == MicrographFamilyState.Available)
				return count;
			count++;
		}
		return -1;
	}

	public void resetFamilyData(MicrographFamilyData mfd)
	{
		String block;
		MetaData emptymd = new MetaData();
		try
		{
			// just in case of user reset
			if (this instanceof SupervisedParticlePicker)
				new File(((SupervisedParticlePicker) this).getTrainingAutoFeaturesVectorFile(mfd)).delete();

			if (!mfd.getAutomaticParticles().isEmpty())
			{
				// removing automatic particles
				block = String.format("%s@%s", mfd.getFamily().getName(), getOutputPath(mfd.getMicrograph().getAutoFilename()));

				emptymd.writeBlock(block);
			}
			if (!mfd.getManualParticles().isEmpty())
			{
				// removing manual particles
				block = String.format("%s@%s", mfd.getFamily().getName(), mfd.getMicrograph().getOFilename());
				emptymd.writeBlock(block);
			}
			mfd.reset();// Resetting family data

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void saveData()
	{
		persistFamilies();
		persistMicrographs();

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
	
	public void exportData(Family f) {
		
		try {
			MetaData md;
			MicrographFamilyData mfd;
			boolean append = false;
			long id;
			String file = getOutputPath(f.getName() + "_extract_list.xmd");
			for (TrainingMicrograph m : micrographs) {

				mfd = m.getFamilyData(f);
				if (!mfd.isEmpty()) {
					md = new MetaData();
					for (TrainingParticle p : mfd.getParticles()) {
						id = md.addObject();
						md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
						md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
						md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					}
					if(!append)
						md.write("mic_" + m.getName() + "@" + file);
					else
						md.writeBlock("mic_" + m.getName() + "@" + file);
					append = true;
				}
			}
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

}
