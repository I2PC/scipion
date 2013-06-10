package xmipp.viewer.particlepicker;

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
import xmipp.jni.PickingClassifier;
import xmipp.utils.TasksManager;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.training.gui.SingleParticlePickerJFrame;
import xmipp.viewer.particlepicker.training.model.AutomaticParticle;
import xmipp.viewer.particlepicker.training.model.MicrographState;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class SingleParticlePicker extends ParticlePicker
{

	protected List<TrainingMicrograph> micrographs;
	private TrainingMicrograph micrograph;
	public static final int defAutoPickPercent = 90;
	private int autopickpercent = defAutoPickPercent;

	private static int mintraining = 15;
	private static String trainingfn = "training.txt";
	private static String trainingmaskfn = "mask.xmp";
	private static String autofeaturesvectorfn = "auto_feature_vectors";
	private static String model = "model";
	private int threads;
	private boolean fastmode;
	private boolean incore;

	public static final int dtemplatesnum = 1;
	private Integer templatesNumber;
	private ImageGeneric templates;
	private String templatesfile;
	private int templateindex;
	private PickingClassifier classifier;

	public SingleParticlePicker(String selfile, String outputdir, Mode mode)
	{
		this(null, selfile, outputdir, mode);

	}

	public SingleParticlePicker(String block, String selfile, String outputdir, Mode mode)
	{

		super(block, selfile, outputdir, mode);
		try
		{
			templatesfile = getOutputPath("templates.stk");
			if (!new File(templatesfile).exists())
			{
				templatesNumber = dtemplatesnum;
				initTemplates();
			}
			else
			{
				this.templates = new ImageGeneric(templatesfile);
				templatesNumber = ((int) templates.getNDim());
				templates.read(templatesfile, false);
			}
			for (TrainingMicrograph m : micrographs)
				loadMicrographData(m);
			System.out.println(getSize());
			classifier = new PickingClassifier(getSize(), getOutputPath("model"));
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new IllegalArgumentException();
		}
	}

	public SingleParticlePicker(String selfile, String outputdir, String reviewfile)
	{
		this(selfile, outputdir, Mode.Review);
		if (!new File(reviewfile).exists())
			throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("review file", reviewfile));
		importAllParticles(reviewfile);
	}

	public SingleParticlePicker(String selfile, String outputdir, Integer threads, boolean fastmode, boolean incore)
	{
		this(selfile, outputdir, Mode.Supervised);

		this.threads = threads;
		this.fastmode = fastmode;
		this.incore = incore;

	}

	public void saveData()
	{
		super.saveData();
		for (Micrograph m : micrographs)
			saveData(m);
		setChanged(false);
	}

	public boolean isFastMode()
	{
		return fastmode;
	}

	public boolean isIncore()
	{
		return incore;
	}

	public int getThreads()
	{
		return threads;
	}

	public synchronized void initTemplates()
	{
		if (templatesNumber == 0)
			return;
		try
		{
			this.templates = new ImageGeneric(ImageGeneric.Float);
			templates.resize(getSize(), getSize(), 1, templatesNumber);
			templates.write(templatesfile);
			templateindex = 0;
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public int getTemplatesNumber()
	{
		if (templatesNumber == null)
			templatesNumber = dtemplatesnum;
		return templatesNumber;
	}

	public void setTemplatesNumber(int num)
	{
		if (num <= 0)
			throw new IllegalArgumentException(
					XmippMessage.getIllegalValueMsgWithInfo("Templates Number", Integer.valueOf(num), "Must have at least one template"));

		this.templatesNumber = num;
		initUpdateTemplates();
	}

	public void setSize(int size)
	{
		super.setSize(size);
		classifier.setSize(size);
		initUpdateTemplates();
	}

	public synchronized ImageGeneric getTemplates()
	{
		return templates;
	}

	public synchronized void setTemplate(ImageGeneric ig)
	{
		float[] matrix;
		try
		{
			//TODO getArrayFloat and setArrayFloat must be call from C both in one function
			matrix = ig.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.FIRST_SLICE);
			templates.setArrayFloat(matrix, ImageGeneric.FIRST_IMAGE + templateindex, ImageGeneric.FIRST_SLICE);
			//			System.out.printf("setTemplate " + templateindex + "\n");
			templateindex++;

		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public String getTemplatesFile()
	{
		return templatesfile;
	}

	public void saveTemplates()
	{
		try
		{
			templates.write(getTemplatesFile());
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

	@Override
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
				micrograph = new TrainingMicrograph(filename, psd, ctf);
				micrographs.add(micrograph);
			}
			if (micrographs.isEmpty())
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", getMicrographsSelFile()));
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	@Override
	public List<TrainingMicrograph> getMicrographs()
	{
		return micrographs;
	}

	public void saveData(Micrograph m)
	{
		TrainingMicrograph tm = (TrainingMicrograph) m;
		long id;
		try
		{
			MetaData md = new MetaData();
			String file;
			file = getOutputPath(m.getPosFile());
			if (!m.hasData())
				new File(file).delete();
			else
			{
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_MICROGRAPH_STATE, tm.getState().toString(), id);
				md.setValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, getAutopickpercent(), id);
				md.writeBlock("header@" + file);
				md.destroy();
				md = new MetaData();
				for (TrainingParticle p : tm.getManualParticles())
				{
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
				}

				md.writeBlock(getParticlesBlock(file));//to append
				md.destroy();

			}
			saveAutomaticParticles(tm);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void saveAutomaticParticles(TrainingMicrograph m)
	{
		try
		{

			long id;
			if (m.hasAutomaticParticles())
			{
				String file = getOutputPath(m.getAutoPosFile());
				MetaData md = new MetaData();
				for (AutomaticParticle p : m.getAutomaticParticles())
				{
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
					md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					md.setValueInt(MDLabel.MDL_ENABLED, (!p.isDeleted()) ? 1 : -1, id);
				}
				md.write(file);
				md.destroy();
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public TrainingMicrograph getMicrograph()
	{
		return micrograph;
	}

	@Override
	public void setMicrograph(Micrograph m)
	{
		if (m == null)
			return;
		if (m.equals(micrograph))
			return;
		this.micrograph = (TrainingMicrograph) m;
	}

	@Override
	public boolean isValidSize(int size)
	{
		for (TrainingParticle p : micrograph.getParticles())
			if (!micrograph.fits(p.getX(), p.getY(), size))
				return false;
		return true;
	}

	public void loadConfig()
	{
		super.loadConfig();
		String file = configfile;
		if (!new File(file).exists())
			return;

		try
		{
			MetaData md = new MetaData(file);
			boolean hasautopercent = md.containsLabel(MDLabel.MDL_PICKING_AUTOPICKPERCENT);
			for (long id : md.findObjects())
			{
				if (hasautopercent)
					autopickpercent = md.getValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, id);
				templatesNumber = md.getValueInt(MDLabel.MDL_PICKING_TEMPLATES, id);
				if (templatesNumber == null || templatesNumber == 0)
					templatesNumber = 1;//for compatibility with previous projects
				mode = Mode.valueOf(md.getValueString(MDLabel.MDL_PICKING_STATE, id));

			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	protected void saveConfig(MetaData md, long id)
	{
		try
		{
			super.saveConfig(md, id);
			md.setValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, getAutopickpercent(), id);

			md.setValueInt(MDLabel.MDL_PICKING_TEMPLATES, getTemplatesNumber(), id);
			md.setValueString(MDLabel.MDL_PICKING_STATE, mode.toString(), id);

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void setAutopickpercent(int autopickpercent)
	{
		this.autopickpercent = autopickpercent;
	}

	public int getAutopickpercent()
	{
		return autopickpercent;
	}

	public void setMode(Mode mode)
	{
		if (mode == Mode.Supervised && getManualParticlesNumber() < mintraining)
			throw new IllegalArgumentException(String.format("You should have at least %s particles to go to %s mode", mintraining, Mode.Supervised));
		this.mode = mode;

	}

	public int getManualParticlesNumber()
	{
		int count = 0;
		for (TrainingMicrograph m : micrographs)
			count += m.getManualParticles().size();
		return count;
	}

	public void exportParticles(String file)
	{

		try
		{
			MetaData md = new MetaData();
			boolean append = false;
			long id;
			for (TrainingMicrograph m : micrographs)
			{
				if (m.isEmpty())
					continue;
				for (TrainingParticle p : m.getParticles())
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
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	public int getAutomaticParticlesNumber(double threshold)
	{
		int count = 0;
		for (TrainingMicrograph m : micrographs)
			count += m.getAutomaticParticlesNumber(threshold);
		return count;
	}

	public void loadAutomaticParticles(TrainingMicrograph m, MetaData md, boolean imported)
	{

		int x, y;
		AutomaticParticle particle;
		Double cost;
		boolean deleted;
		try
		{

			for (long id : md.findObjects())
			{

				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				cost = md.getValueDouble(MDLabel.MDL_COST, id);
				if (cost == null)
					throw new IllegalArgumentException("Invalid format for " + md.getFilename());
				deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? false : true;
				particle = new AutomaticParticle(x, y, this, m, cost, deleted);
				m.addAutomaticParticle(particle, imported);
			}
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void loadAutomaticParticles(TrainingMicrograph m, String file, boolean imported)
	{
		if (!new File(file).exists())
			return;
		MetaData md = new MetaData(file);
		loadAutomaticParticles(m, md, imported);
		md.destroy();
	}

	//returns micrograph particles filepath, according to format
	public String getImportMicrographName(String path, String filename, Format f)
	{

		String base = Filename.removeExtension(Filename.getBaseName(filename));
		switch (f)
		{
		case Xmipp24:
			return Filename.join(path, base, base + ".raw.Common.pos");
		case Xmipp30:
			return Filename.join(path, base + ".pos");
		case Xmipp301:
			return Filename.join(path, base + ".pos");
		case Eman:
			return Filename.join(path, base + ".box");

		default:
			return null;
		}

	}

	/** Return the number of particles imported */
	public String importParticlesFromFolder(String path, Format f, float scale, boolean invertx, boolean inverty)
	{
		if (f == Format.Auto)
			f = detectFormat(path);
		if (f == Format.Unknown)
			throw new IllegalArgumentException(XmippMessage.getIllegalValueMsg("format", Format.Unknown));

		String filename;
		String result = "";
		initTemplates();
		for (TrainingMicrograph m : micrographs)
		{
			filename = getImportMicrographName(path, m.getFile(), f);
			//			System.out.println("  filename: " + filename);
			if (Filename.exists(filename))
				result += importParticlesFromFile(filename, f, m, scale, invertx, inverty);
		}

		return result;
	}// function importParticlesFromFolder

	/** Return the number of particles imported from a file */
	public String importParticlesFromFile(String path, Format f, Micrograph m, float scale, boolean invertx, boolean inverty)
	{
		MetaData md = new MetaData();
		fillParticlesMdFromFile(path, f, m, md, scale, invertx, inverty);
		String result = (md != null) ? importParticlesFromMd(m, md) : "";
		saveData(m);
		md.destroy();
		return result;
	}// function importParticlesFromFile

	public String importParticlesFromMd(Micrograph m, MetaData md)
	{

		m.reset();
		TrainingMicrograph tm = (TrainingMicrograph) m;
		long[] ids = md.findObjects();
		int x, y;
		double cost;
		boolean hasCost = md.containsLabel(MDLabel.MDL_COST);
		String result = "";
		int size = getSize();

		for (long id : ids)
		{
			x = md.getValueInt(MDLabel.MDL_XCOOR, id);
			y = md.getValueInt(MDLabel.MDL_YCOOR, id);
			if (!m.fits(x, y, size))// ignore out of bounds particle
			{
				result += XmippMessage.getOutOfBoundsMsg("Particle") + String.format(" on x:%s y:%s", x, y);
				continue;
			}
			cost = hasCost ? md.getValueDouble(MDLabel.MDL_COST, id) : 0;
			if (cost == 0 || cost > 1)
				tm.addManualParticle(new TrainingParticle(x, y, this, tm, cost), this, false, true);
			else
				tm.addAutomaticParticle(new AutomaticParticle(x, y, this, tm, cost, false), true);
		}
		return result;
	}// function importParticlesFromMd

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
			initTemplates();
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
			saveData();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
		return null;

	}// function importAllParticles

	public synchronized void resetParticleImages()//to update templates with the right particles
	{
		for (TrainingMicrograph m : micrographs)
		{
			for (TrainingParticle p : m.getManualParticles())
				p.resetImagePlus();

		}
	}

	public void initUpdateTemplates()
	{
		TasksManager.getInstance().addTask(new UpdateTemplatesTask(this));
		saveConfig();
	}

	public synchronized void updateTemplates()
	{
		try
		{
			if (getMode() != Mode.Manual)
				return;

			initTemplates();
			ImageGeneric igp;
			List<TrainingParticle> particles;
			TrainingParticle particle;
			double[] align;

			for (TrainingMicrograph m : getMicrographs())
			{
				for (int i = 0; i < m.getManualParticles().size(); i++)
				{
					particles = m.getManualParticles();
					particle = particles.get(i);
					igp = particle.getImageGeneric();
					if (getTemplateIndex() < getTemplatesNumber())
						setTemplate(igp);
					else
					{
						align = getTemplates().alignImage(igp);
						applyAlignment(particle, igp, align);
					}
				}
			}
			saveTemplates();
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public synchronized void addParticleToTemplates(TrainingParticle particle)
	{
		try
		{
			ImageGeneric igp = particle.getImageGeneric();
			// will happen only in manual mode
			if (getTemplateIndex() < getTemplatesNumber())
				setTemplate(igp);
			else
			{
				double[] align = getTemplates().alignImage(igp);
				applyAlignment(particle, igp, align);
			}
			saveTemplates();
			//			System.out.println("templates saved");

		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}

	}

	public boolean hasManualParticles()
	{
		for (TrainingMicrograph m : micrographs)
		{
			if (m.hasManualParticles())
				return true;
		}
		return false;
	}

	public int getTemplateIndex()
	{
		return templateindex;
	}

	public synchronized void centerParticle(TrainingParticle p)
	{
		if (getManualParticlesNumber() <= templatesNumber)
			return;//templates not ready
		Particle shift = null;
		try
		{
			ImageGeneric igp = p.getImageGeneric();
			shift = templates.bestShift(igp);
			double distance = Math.sqrt(Math.pow(shift.getX(), 2) + Math.pow(shift.getY(), 2)) / getSize();
			System.out.printf("normalized distance:%.2f\n", distance);
			if (distance < 0.25)
			{
				p.setX(p.getX() + shift.getX());
				p.setY(p.getY() + shift.getY());
			}

		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public synchronized void applyAlignment(TrainingParticle particle, ImageGeneric igp, double[] align)
	{
		try
		{
			particle.setLastalign(align);
			templates.applyAlignment(igp, particle.getTemplateIndex(), particle.getTemplateRotation(), particle.getTemplateTilt(), particle
					.getTemplatePsi());
			//System.out.printf("adding particle: %d %.2f %.2f %.2f\n", particle.getTemplateIndex(), particle.getTemplateRotation(), particle.getTemplateTilt(), particle.getTemplatePsi());
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void loadMicrographData(TrainingMicrograph micrograph)
	{
		try
		{
			String file = getOutputPath(micrograph.getPosFile());
			MicrographState state;
			Integer autopickpercent;
			if (!new File(file).exists())
				return;

			MetaData md = new MetaData("header@" + file);
			boolean hasautopercent = md.containsLabel(MDLabel.MDL_PICKING_AUTOPICKPERCENT);
			for (long id : md.findObjects())
			{
				state = MicrographState.valueOf(md.getValueString(MDLabel.MDL_PICKING_MICROGRAPH_STATE, id));
				if (hasautopercent)
					autopickpercent = md.getValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, id);
				else
					autopickpercent = 50;// compatibility with previous projects
				micrograph.setState(state);
				micrograph.setAutopickpercent(autopickpercent);

			}
			md.destroy();
			loadManualParticles(micrograph, file);
			loadAutomaticParticles(micrograph, getOutputPath(micrograph.getAutoPosFile()), false);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void loadManualParticles(TrainingMicrograph micrograph, String file)
	{
		if (!new File(file).exists())
			return;

		int x, y;
		TrainingParticle particle;

		try
		{
			MetaData md = new MetaData(getParticlesBlock(file));

			for (long id : md.findObjects())
			{
				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				particle = new TrainingParticle(x, y, this, micrograph);
				micrograph.addManualParticle(particle, this, false, false);
			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void train(SingleParticlePickerJFrame frame)
	{
		if (mode != Mode.Manual)
			throw new IllegalArgumentException(XmippMessage.getIllegalStateForOperationMsg("picker", this.mode.toString()));
		frame.getCanvas().setEnabled(false);
		XmippWindowUtil.blockGUI(frame, "Training...");
		MetaData trainmd = new MetaData();
		MetaData outputmd = new MetaData();

		micrograph.setState(MicrographState.Supervised);
		saveData();
		long id;
		String posfile;
		for (TrainingMicrograph m : micrographs)
		{
			if (m.hasManualParticles())
			{
				posfile = getOutputPath(m.getPosFile());
				id = trainmd.addObject();
				trainmd.setValueString(MDLabel.MDL_MICROGRAPH, m.getFile(), id);
				trainmd.setValueString(MDLabel.MDL_MICROGRAPH_PARTICLES, posfile, id);
				m.getAutomaticParticles().clear();
				new File(getOutputPath(m.getAutoPosFile())).delete();

			}
		}


		new Thread(new TrainRunnable(frame, trainmd, outputmd)).start();
//		new TrainRunnable(frame, trainmd, outputmd).run();
	}

	public class TrainRunnable implements Runnable
	{

		private SingleParticlePickerJFrame frame;
		private MetaData trainmd;
		private MetaData outputmd;

		public TrainRunnable(SingleParticlePickerJFrame frame, MetaData trainmd, MetaData outputmd)
		{
			this.frame = frame;
			this.trainmd = trainmd;
			this.outputmd = outputmd;

		}

		public void run()
		{
			try
			{
				System.out.println("Train started");
				classifier.train(trainmd);//should remove training files
				micrograph.setAutopickpercent(autopickpercent);
				classifier.autopick(micrograph.getFile(), outputmd, micrograph.getAutopickpercent());
				loadAutomaticParticles(micrograph, outputmd, false);
				String path = getParticlesBlock(getOutputPath(micrograph.getAutoPosFile()));
				outputmd.write(path);
				frame.getCanvas().repaint();
				XmippWindowUtil.releaseGUI(frame.getRootPane());
				frame.getCanvas().setEnabled(true);
				frame.updateMicrographsModel();
				trainmd.destroy();
				outputmd.destroy();
			}
			catch (Exception e)
			{
				ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
				throw new IllegalArgumentException(e.getMessage());
			}
		}
	}

	public void autopick(SingleParticlePickerJFrame frame)
	{
		micrograph.setState(MicrographState.Supervised);
		saveData(micrograph);

		MetaData outputmd = new MetaData();
		frame.getCanvas().setEnabled(false);
		XmippWindowUtil.blockGUI(frame, "Autopicking...");

		new Thread(new AutopickRunnable(frame, outputmd)).start();

	}

	public class AutopickRunnable implements Runnable
	{

		private SingleParticlePickerJFrame frame;
		private MetaData outputmd;

		public AutopickRunnable(SingleParticlePickerJFrame frame, MetaData outputmd)
		{
			this.frame = frame;
			this.outputmd = outputmd;
		}

		public void run()
		{
			micrograph.getAutomaticParticles().clear();
			micrograph.setAutopickpercent(autopickpercent);
			classifier.autopick(micrograph.getFile(), outputmd, micrograph.getAutopickpercent());
			loadAutomaticParticles(micrograph, outputmd, false);
			String path = getParticlesBlock(getOutputPath(micrograph.getAutoPosFile()));
			outputmd.write(path);
			frame.getCanvas().repaint();
			frame.getCanvas().setEnabled(true);
			frame.updateMicrographsModel();
			XmippWindowUtil.releaseGUI(frame.getRootPane());
			outputmd.destroy();
		}

	}

	public void correct()
	{
		getMicrograph().setState(MicrographState.Corrected);
		saveData(micrograph);
		MetaData manualmd = new MetaData(getOutputPath(micrograph.getPosFile()));
		MetaData automaticmd = new MetaData(getOutputPath(micrograph.getAutoPosFile()));

		new Thread(new CorrectRunnable(manualmd, automaticmd)).start();

	}

	public class CorrectRunnable implements Runnable
	{

		private MetaData manualmd;
		private MetaData automaticmd;

		public CorrectRunnable(MetaData manualmd, MetaData automaticmd)
		{
			this.manualmd = manualmd;
			this.automaticmd = automaticmd;
		}

		public void run()
		{

			classifier.correct(manualmd, automaticmd);
			manualmd.destroy();
			automaticmd.destroy();
		}
	}

}
