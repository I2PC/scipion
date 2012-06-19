package xmipp.particlepicker.training.model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;

import xmipp.particlepicker.Family;
import xmipp.particlepicker.ParticlePicker;
import xmipp.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.utils.XmippMessage;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

public abstract class TrainingPicker extends ParticlePicker
{

	protected List<TrainingMicrograph> micrographs;

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
		this.micrographs = new ArrayList<TrainingMicrograph>();

	}
	
	public TrainingPicker(String selfile, String outputdir, FamilyState mode)
	{
		super(selfile, outputdir, mode);
		this.micrographs = new ArrayList<TrainingMicrograph>();

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
				if (md.containsLabel(MDLabel.MDL_MICROGRAPH))
					file = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
				else if (md.containsLabel(MDLabel.MDL_IMAGE))
					file = md.getValueString(MDLabel.MDL_IMAGE, id);
				else
					throw new IllegalArgumentException(String.format("Labels MDL_MICROGRAPH or MDL_IMAGE not found in metadata %s", selfile));

				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new TrainingMicrograph(file, ctf, families, getMode());
				loadMicrographData(micrograph);
				
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", selfile));
		

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
			String file = getOutputPath(micrograph.getPosFile());
			String fname;
			Family family;
			MicrographFamilyState state;
			MicrographFamilyData mfd;
			List<MicrographFamilyData> mfdatas = new ArrayList<MicrographFamilyData>();
			if (!new File(file).exists())
				return;
			MetaData md = new MetaData("families@" + file);
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
				loadManualParticles(mfd, file);
				file = getOutputPath(micrograph.getAutoPosFile());
				loadAutomaticParticles(mfd, file, false);
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
		
		MetaData md;
		try
		{
			md = new MetaData(family.getName() + "@" + file);

			for (long id : md.findObjects())
			{

				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				particle = new TrainingParticle(x, y, family, mfd.getMicrograph());
				mfd.addManualParticle(particle);
			}
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
		MetaData md;
		boolean deleted;
		try
		{
			md = new MetaData(f.getName() + "@" + file);

			for (long id : md.findObjects())
			{

				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				cost = md.getValueDouble(MDLabel.MDL_COST, id);
				if(cost == null)
					throw new IllegalArgumentException("Invalid format for " + file);
				System.out.println(cost);
				deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? false : true;
				particle = new AutomaticParticle(x, y, f, mfd.getMicrograph(), cost, deleted);
				mfd.addAutomaticParticle(particle, imported);
			}
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
			String file;
			for (TrainingMicrograph m : micrographs)
			{
				file = getOutputPath(m.getPosFile());
				if (!m.hasData())
					new File(file).delete();
				else
				{
					persistMicrographFamilies(m);
					for (MicrographFamilyData mfd : m.getFamiliesData())
					{
						md = new MetaData();
						for (TrainingParticle p : mfd.getManualParticles())
						{
							id = md.addObject();
							md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
							md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
						}
						block = mfd.getFamily().getName() + "@" + file;
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
			new File(getOutputPath(m.getAutoPosFile())).delete();
		else
			for (MicrographFamilyData mfd : m.getFamiliesData())
				persistAutomaticParticles(mfd);
	}

	
	public void persistAutomaticParticles(MicrographFamilyData mfd)
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
					md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
					md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					md.setValueInt(MDLabel.MDL_ENABLED, (!p.isDeleted()) ? 1 : -1, id);
				}
				md.write(section);
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
			String file = getOutputPath(m.getPosFile());
			MetaData md = new MetaData();
			for (MicrographFamilyData mfd : m.getFamiliesData())
			{
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, mfd.getFamily().getName(), id);
				md.setValueString(MDLabel.MDL_PICKING_MICROGRAPH_FAMILY_STATE, mfd.getState().toString(), id);
			}
			md.writeBlock("families@" + file);

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
			persistMicrographs();
		}
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
			MetaData md;
			MicrographFamilyData mfd;
			boolean append = false;
			long id;
			for (TrainingMicrograph m : micrographs)
			{

				mfd = m.getFamilyData(family);
				if (!mfd.isEmpty())
				{
					md = new MetaData();
					for (TrainingParticle p : mfd.getParticles())
					{
						id = md.addObject();
						md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
						md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
						md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					}
					if (!append)
						md.write("mic_" + m.getName() + "@" + file);
					else
						md.writeBlock("mic_" + m.getName() + "@" + file);
					append = true;
				}
			}
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	@Override
	public void importParticlesFromXmipp24Folder(String path)
	{
		throw new UnsupportedOperationException(XmippMessage.getNotImplementedYetMsg());

	}
	
	@Override
	public void importParticlesFromXmipp30Folder(String dir)
	{
		MicrographFamilyData mfd;
		for(TrainingMicrograph m: micrographs)
		{
			mfd = m.getFamilyData(family);
			loadManualParticles(mfd, dir + File.separator + m.getPosFile());
			loadAutomaticParticles(mfd, dir + File.separator + m.getAutoPosFile(), true);//boolean for imported, so that available micrographs can have automatic particles
		}
	}

	
	
	public void importParticlesFromEmanFolder(String folder)
	{
		String file;
		for (TrainingMicrograph tm : micrographs)
		{
			file = folder + "/" + tm.getName() + ".box";
			importParticlesFromEmanFile(tm.getFamilyData(family), file);
		}
	}

	public void importParticlesFromXmipp24File(MicrographFamilyData familyData, String file)
	{
		throw new UnsupportedOperationException(XmippMessage.getNotImplementedYetMsg());

	}

	public void importParticlesFromXmipp30File(MicrographFamilyData mfd, String file)
	{
		try
		{
			loadAutomaticParticles(mfd, file, true);//if manual raises exception
		}
		catch(Exception e)
		{
			loadManualParticles(mfd, file);
		}
	}
	
	
	public void importParticlesFromEmanFile(MicrographFamilyData mfd, String file)
	{
		if(!new File(file).exists())
			return;
		BufferedReader reader = null;
		try
		{
			reader = new BufferedReader(new FileReader(file));
			boolean inverty = false;
			String line = reader.readLine();
			reader.close();
			if(line.split("\t").length > 4)//eman 1.0
				inverty = true;
			MetaData md = new MetaData();
			md.readPlain(file, "Xcoor Ycoor particleSize");
			long[] ids;
			int x, y, size = 0, height;
			Double cost = 2.0;
			ids = md.findObjects();
			for (long id : ids)
			{
				size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
				x = md.getValueInt(MDLabel.MDL_XINT, id) + size / 2;
				y = md.getValueInt(MDLabel.MDL_YINT, id) + size / 2;
				if(inverty)
				{
					height = mfd.getMicrograph().getImagePlus().getHeight();
					y = height - y;
				}
				mfd.addManualParticle(new TrainingParticle(x, y, mfd.getFamily(), mfd.getMicrograph(), cost));

			}
			if (size > 0)
				mfd.getFamily().setSize(size);

		}
		catch (Exception e)
		{
			if(reader != null)
				try {
					reader.close();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}
	
	public void importAllParticles(String file)
	{// Expected a file for all micrographs
		try
		{
			MetaData md;
			long[] ids;
			int x, y;
			Double cost;

			List<String> blocks = Arrays.asList(MetaData.getBlocksInMetaDataFile(file));
			String block;
			for (TrainingMicrograph m : micrographs)
			{
				m.reset();

				block = "mic_" + m.getName();
				if (blocks.contains(block))
				{
					md = new MetaData(block + "@" + file);

					ids = md.findObjects();
					for (long id : ids)
					{

						x = md.getValueInt(MDLabel.MDL_XINT, id);
						y = md.getValueInt(MDLabel.MDL_YINT, id);
						cost = md.getValueDouble(MDLabel.MDL_COST, id);
						if (cost == null || cost == 0 || cost > 1)
							m.addManualParticle(new TrainingParticle(x, y, family, m, cost));
						else
							m.addAutomaticParticle(new AutomaticParticle(x, y, family, m, cost, false), true);
					}
				}
			}
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void removeFamily(Family family)
	{
		if (getManualParticlesNumber(family) > 0)// perhaps I have to check
													// automatic particles
			throw new IllegalArgumentException(XmippMessage.getAssociatedDataMsg("family"));
		if (families.size() == 1)
			throw new IllegalArgumentException(XmippMessage.getIllegalDeleteMsg("family"));
		families.remove(family);
		for(TrainingMicrograph m: micrographs)
			m.removeFamilyData(family);
	}
	


}
