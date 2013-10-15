package xmipp.viewer.particlepicker.training.model;

import java.awt.Point;
import java.awt.Rectangle;
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
import xmipp.utils.DEBUG;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.training.gui.SingleParticlePickerJFrame;

/**
 * Business object for Single Particle Picker GUI. Inherits from ParticlePicker
 * 
 * @author airen
 * 
 */
public class SingleParticlePicker extends ParticlePicker {

	protected List<SingleParticlePickerMicrograph> micrographs;
	private SingleParticlePickerMicrograph micrograph;
	public static final int defAutoPickPercent = 90;
	private int autopickpercent = defAutoPickPercent;

	public static final int mintraining = 15;

	private int threads = 1;
	private boolean fastmode = true;
	private boolean incore = false;

	public static int dtemplatesnum = 1;
	private ImageGeneric templates;
	private ImageGeneric radialtemplates;
	private String templatesfile;
	private int templateindex;
	private PickingClassifier classifier;

	// private String reviewfile;

	public SingleParticlePicker(String selfile, String outputdir, Mode mode) {
		this(null, selfile, outputdir, mode);

	}

	public SingleParticlePicker(String block, String selfile, String outputdir,
			Mode mode) {

		super(block, selfile, outputdir, mode);
		try {
			templatesfile = getOutputPath("templates.stk");
			if (!new File(templatesfile).exists())

				initTemplates(dtemplatesnum);
			else {
				this.templates = new ImageGeneric(templatesfile);
				templates.read(templatesfile, false);
				radialtemplates = new ImageGeneric(ImageGeneric.Float);
				radialtemplates.resize(getSize(), getSize(), 1,
						getTemplatesNumber());
			}

			templates.getRadialAvg(radialtemplates);

			for (SingleParticlePickerMicrograph m : micrographs)
				loadMicrographData(m);
			classifier = new PickingClassifier(getSize(),
					getOutputPath("model"));
		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException();
		}
	}

	// public SingleParticlePicker(String selfile, String outputdir,
	// String reviewfile) {
	// this(selfile, outputdir, Mode.Review);
	// if (!new File(reviewfile).exists())
	// throw new IllegalArgumentException(
	// XmippMessage.getNoSuchFieldValueMsg("review file",
	// reviewfile));
	// this.reviewfile = reviewfile;
	// importAllParticles(reviewfile);
	// }

	public SingleParticlePicker(String selfile, String outputdir,
			Integer threads, boolean fastmode, boolean incore) {
		this(selfile, outputdir, Mode.Manual);

		this.threads = threads;
		this.fastmode = fastmode;
		this.incore = incore;

	}

	/* Save changes to the current micrograph. */
	public void saveData() {
		super.saveData();
		saveData(micrograph);
		setChanged(false);
	}

	/* Save changes to ALL micrographs, mainly used when importing. */
	public void saveAllData() {
		super.saveData();
		for (Micrograph m : micrographs)
			saveData(m);
		setChanged(false);
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

	public synchronized void initTemplates() {
		initTemplates(getTemplatesNumber());
	}

	public synchronized void initTemplates(int num) {
		System.out.println("initTemplates");
		if (num == 0)
			return;
		try {
			this.templates = new ImageGeneric(ImageGeneric.Float);
			templates.resize(getSize(), getSize(), 1, num);
			templates.write(templatesfile);
			templateindex = 0;
			radialtemplates = new ImageGeneric(ImageGeneric.Float);
			radialtemplates.resize(getSize(), getSize(), 1, num);
		} catch (Exception e) {
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public int getTemplatesNumber() {
		if (templates == null)
			return dtemplatesnum;
		try {
			return (int) templates.getNDim();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return dtemplatesnum;
	}

	public synchronized void setTemplatesNumber(int num) {
		System.out.println("setTemplatesNumber");
		if (num <= 0)
			throw new IllegalArgumentException(
					XmippMessage.getIllegalValueMsgWithInfo("Templates Number",
							Integer.valueOf(num),
							"Must have at least one template"));

		initUpdateTemplates(num);
	}

	public void setSize(int size) {

		super.setSize(size);
		classifier.setSize(size);
	}

	public synchronized ImageGeneric getTemplates() {
		System.out.println("getTemplates");
		return templates;
	}

	public synchronized void setTemplate(ImageGeneric ig) {
		System.out.println("setTemplate");
		float[] matrix;
		try {
			// TODO getArrayFloat and setArrayFloat must be call from C both in
			// one function
			matrix = ig.getArrayFloat(ImageGeneric.FIRST_IMAGE,
					ImageGeneric.FIRST_SLICE);
			templates.setArrayFloat(matrix, ImageGeneric.FIRST_IMAGE
					+ templateindex, ImageGeneric.FIRST_SLICE);
			// System.out.printf("setTemplate " + templateindex + "\n");
			templateindex++;

		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public String getTemplatesFile() {
		return templatesfile;
	}

	public void saveTemplates() {
		try {
			templates.write(getTemplatesFile());
		} catch (Exception e) {
			throw new IllegalArgumentException(e);
		}

	}

	@Override
	public void loadEmptyMicrographs() {
		if (micrographs == null)
			micrographs = new ArrayList<SingleParticlePickerMicrograph>();
		else
			micrographs.clear();
		SingleParticlePickerMicrograph micrograph;
		String psd = null, ctf = null, filename;
		try {
			MetaData md = new MetaData(getMicrographsSelFile());
			md.removeDisabled();
			int fileLabel;

			if (md.containsLabel(MDLabel.MDL_MICROGRAPH))
				fileLabel = MDLabel.MDL_MICROGRAPH;
			else if (md.containsLabel(MDLabel.MDL_IMAGE))
				fileLabel = MDLabel.MDL_IMAGE;
			else
				throw new IllegalArgumentException(
						String.format(
								"Labels MDL_MICROGRAPH or MDL_IMAGE not found in metadata %s",
								selfile));
			boolean existspsd = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			boolean existsctf = md.containsLabel(MDLabel.MDL_CTF_MODEL);
			long[] ids = md.findObjects();
			for (long id : ids) {

				filename = md.getValueString(fileLabel, id);
				if (existspsd)
					psd = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_CTF_MODEL, id);
				micrograph = new SingleParticlePickerMicrograph(filename, psd,
						ctf);
				micrographs.add(micrograph);
			}
			if (micrographs.isEmpty())
				throw new IllegalArgumentException(String.format(
						"No micrographs specified on %s",
						getMicrographsSelFile()));
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	@Override
	public List<SingleParticlePickerMicrograph> getMicrographs() {
		return micrographs;
	}

	public void saveData(Micrograph m) {

		SingleParticlePickerMicrograph tm = (SingleParticlePickerMicrograph) m;
		long id;
		try {
			MetaData md = new MetaData();
			String file;
			file = getOutputPath(m.getPosFile());
			if (!m.hasData())
				new File(file).delete();
			else {
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_MICROGRAPH_STATE, tm
						.getState().toString(), id);
				md.setValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT,
						getAutopickpercent(), id);
				md.setValueDouble(MDLabel.MDL_COST, tm.getThreshold(), id);
				md.write("header@" + file); // using md.write will remove other
											// existing blocks
				md.clear();
				
				for (ManualParticle p : tm.getManualParticles()) {
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);

				}
				md.writeBlock(getParticlesBlock(file));// to replace/append
			}
			md.destroy();
			saveAutomaticParticles(tm);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void saveAutomaticParticles(SingleParticlePickerMicrograph m) {
		try {

			long id;
			if (m.hasAutomaticParticles()) {
				MetaData md = new MetaData();
				for (AutomaticParticle p : m.getAutomaticParticles()) {
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
					md.setValueDouble(MDLabel.MDL_COST, p.getCost(), id);
					md.setValueInt(MDLabel.MDL_ENABLED, (!p.isDeleted()) ? 1
							: -1, id);
				}
				md.writeBlock(getParticlesAutoBlock(m));// to replace/append
				md.destroy();
			}

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public SingleParticlePickerMicrograph getMicrograph() {
		return micrograph;
	}

	@Override
	public void setMicrograph(Micrograph m) {
		if (m == null)
			return;
		if (m.equals(micrograph))
			return;
		this.micrograph = (SingleParticlePickerMicrograph) m;
	}

	@Override
	public boolean isValidSize(int size) {
		for (SingleParticlePickerMicrograph m : micrographs)
			for (ManualParticle p : m.getParticles())
				if (!micrograph.fits(p.getX(), p.getY(), size))
					return false;
		return true;
	}

	public void loadConfig() {
		super.loadConfig();
		String file = configfile;
		if (!new File(file).exists())
			return;

		try {
			MetaData md = new MetaData(file);
			Mode configmode;
			boolean hasautopercent = md
					.containsLabel(MDLabel.MDL_PICKING_AUTOPICKPERCENT);
			for (long id : md.findObjects()) {
				if (hasautopercent)
					autopickpercent = md.getValueInt(
							MDLabel.MDL_PICKING_AUTOPICKPERCENT, id);
				dtemplatesnum = md.getValueInt(MDLabel.MDL_PICKING_TEMPLATES,
						id);
				if (dtemplatesnum == 0)
					dtemplatesnum = 1;// for compatibility with previous
										// projects
				configmode = Mode.valueOf(md.getValueString(
						MDLabel.MDL_PICKING_STATE, id));
				if (mode == Mode.Supervised && configmode == Mode.Manual)
					throw new IllegalArgumentException(
							"Cannot switch to Supervised mode from the command line");
				if (mode == Mode.Manual && configmode == Mode.Supervised)
					mode = Mode.Supervised;

			}
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	protected void saveConfig(MetaData md, long id) {
		try {
			super.saveConfig(md, id);
			md.setValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT,
					getAutopickpercent(), id);
			md.setValueInt(MDLabel.MDL_PICKING_TEMPLATES, getTemplatesNumber(),
					id);
			md.setValueString(MDLabel.MDL_PICKING_STATE, mode.toString(), id);

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void setAutopickpercent(int autopickpercent) {
		this.autopickpercent = autopickpercent;
	}

	public int getAutopickpercent() {
		return autopickpercent;
	}

	public void setMode(Mode mode) {
		if (mode == Mode.Supervised && getManualParticlesNumber() < mintraining)
			throw new IllegalArgumentException(String.format(
					"You should have at least %s particles to go to %s mode",
					mintraining, Mode.Supervised));
		this.mode = mode;
		if (mode == Mode.Manual)
			convertAutomaticToManual();

	}

	protected void convertAutomaticToManual() {
		ManualParticle p;
		for (SingleParticlePickerMicrograph m : micrographs) {
			m.deleteBelowThreshold();
			for (AutomaticParticle ap : m.getAutomaticParticles()) {
				if (!ap.isDeleted()) {
					p = new ManualParticle(ap.getX(), ap.getY(), this, m);
					m.addManualParticle(p, this, false);
				}
			}

			m.getAutomaticParticles().clear();
			if (m.hasManualParticles())
				m.setState(MicrographState.Manual);
			else
				m.setState(MicrographState.Available);
		}

		saveAllData();
		initUpdateTemplates();
	}

	public int getManualParticlesNumber() {
		int count = 0;
		for (SingleParticlePickerMicrograph m : micrographs)
			count += m.getManualParticles().size();
		return count;
	}

	public void exportParticles(String file) {

		try {
			MetaData md = new MetaData();
			boolean append = false;
			long id;
			for (SingleParticlePickerMicrograph m : micrographs) {
				if (m.isEmpty())
					continue;
				for (ManualParticle p : m.getParticles()) {
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
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	public int getAutomaticParticlesNumber(double threshold) {
		int count = 0;
		for (SingleParticlePickerMicrograph m : micrographs)
			count += m.getAutomaticParticlesNumber(threshold);
		return count;
	}

	public void loadAutomaticParticles(SingleParticlePickerMicrograph m,
			MetaData md) {

		int x, y;
		AutomaticParticle particle;
		Double cost;
		boolean deleted;
		try {
			boolean added = false;

			for (long id : md.findObjects()) {

				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				cost = md.getValueDouble(MDLabel.MDL_COST, id);
				if (cost == null)
					throw new IllegalArgumentException("Invalid format for "
							+ md.getFilename());
				deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? false
						: true;
				particle = new AutomaticParticle(x, y, this, m, cost, deleted);
				m.addAutomaticParticle(particle);
				added = true;
			}
			// Change the state of the micrograph if some automatic particles
			// were added
			if (added && m.getState() == MicrographState.Available)
				m.setState(MicrographState.Auto);

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void loadAutomaticParticles(SingleParticlePickerMicrograph m,
			String file) {
		if (new File(file).exists()
				&& MetaData.containsBlock(file, particlesAutoBlock))// todo:
																	// check the
																	// auto
																	// block
																	// exists
		{
			MetaData md = new MetaData(getParticlesAutoBlock(file));
			loadAutomaticParticles(m, md);
			md.destroy();
		}
	}

	// returns micrograph particles filepath, according to format
	public String getImportMicrographName(String path, String filename, Format f) {

		String base = Filename.removeExtension(Filename.getBaseName(filename));
		switch (f) {
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
	public String importParticlesFromFolder(String path, Format f, float scale,
			boolean invertx, boolean inverty) {

		//System.err.format("JM_DEBUG:   >>>> importParticlesFromFolder   scale: %f\n", scale);
		if (f == Format.Auto)
			f = detectFormat(path);
		if (f == Format.Unknown)
			throw new IllegalArgumentException("Unable to detect format");

		String filename;
		String result = "";

		String particlesfile = getExportFile(path);
		if (particlesfile != null && f == Format.Xmipp30) {
			importAllParticles(particlesfile);
			return "";
		}
		for (SingleParticlePickerMicrograph m : micrographs) {
			filename = getImportMicrographName(path, m.getFile(), f);
			if (Filename.exists(filename)) {
				result += importParticlesFromFile(filename, f, m, scale,
						invertx, inverty);
			}
		}

		return result;
	}// function importParticlesFromFolder

	/** Return the number of particles imported from a file */
	public String importParticlesFromFile(String path, Format f,
			SingleParticlePickerMicrograph m, float scale, boolean invertx,
			boolean inverty) {
		String result = "";
		try {
			MetaData md = new MetaData();
			fillParticlesMdFromFile(path, f, m, md, scale, invertx, inverty);
			System.err.format("After fillParticlesMdFromFile, size:\n", md.size());
			if (md.size() > 0) {
				//Be sure that width and height are loaded
				result = importParticlesFromMd(m, md);
			}
			md.destroy();
		} catch (Exception ex) {
			System.err.format("Error importing from file '%s' message: '%s'\n",
					path, ex.toString());
		}
		return result;
	}// function importParticlesFromFile

	public String importParticlesFromMd(SingleParticlePickerMicrograph m,
			MetaData md) {

		resetMicrograph(m);
		m.loadDimensions();
		long[] ids = md.findObjects();
		int x, y;
		String result = "";
		int size = getSize();

		for (long id : ids) {
			x = md.getValueInt(MDLabel.MDL_XCOOR, id);
			y = md.getValueInt(MDLabel.MDL_YCOOR, id);
			if (!m.fits(x, y, size))// ignore out of bounds particles
			{
				result += XmippMessage.getOutOfBoundsMsg(String.format(
						"Particle on x:%s y:%s\n", x, y)) + " dismissed\n";

				continue;
			}
			m.addManualParticle(new ManualParticle(x, y, this, m, 2), this,
					false);

		}
		saveData(m);
		return result;
	}// function importParticlesFromMd

	public boolean isReviewFile(String file) {
		if(!file.endsWith("xmd"))
			return false;
		boolean result = false;
		boolean nullMicrograph = false;
		try {
			MetaData md = new MetaData(file);
			
			if (md.containsLabel(MDLabel.MDL_XCOOR) &&
				md.containsLabel(MDLabel.MDL_YCOOR) &&
				md.containsLabel(MDLabel.MDL_COST)) {
				String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
				List<String> blocks = Arrays.asList(blocksArray);
				
				String name;
				for (String block : blocks) {
					name = block.replace("mic_", "");
					if (getMicrograph(name) == null) {
						nullMicrograph = true;
						break;
					}
				}
			}
			md.destroy();

			if (!nullMicrograph)
				result = true;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return result;
	}

	public void importAllParticles(String file) {// Expected a file for all
													// micrographs
		try {
			String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
			List<String> blocks = Arrays.asList(blocksArray);
			String block;
			MetaData md = new MetaData();

			for (SingleParticlePickerMicrograph m : micrographs) {
				resetMicrograph(m);
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

	public String importAllParticles(String file, float scale, boolean invertx,
			boolean inverty) {// Expected a file for all
								// micrographs
		try {

			String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
			List<String> blocks = Arrays.asList(blocksArray);
			String block;
			MetaData md = new MetaData();
			int width, height;

			for (SingleParticlePickerMicrograph m : micrographs) {
				resetMicrograph(m);
				block = "mic_" + m.getName();
				if (blocks.contains(block)) {
					String blockName = block + "@" + file;
					md.read(blockName);
					width = (int) (m.width / scale);// original width
					height = (int) (m.height / scale);// original height
					if (invertx)
						md.operate(String.format("xcoor=%d-xcoor", width));
					if (inverty)
						md.operate(String.format("ycoor=%d-ycoor", height));
					if (scale != 1.f)
						md.operate(String.format(
								"xcoor=xcoor*%f,ycoor=ycoor*%f", scale, scale));
					importParticlesFromMd(m, md);
				}
			}
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
		return null;

	}// function importAllParticles

	public synchronized void resetParticleImages()// to update templates with
													// the right particles
	{
		System.out.println("resetParticleImages");
		for (SingleParticlePickerMicrograph m : micrographs) {
			for (ManualParticle p : m.getManualParticles())
				p.resetImagePlus();

		}
	}

	public void initUpdateTemplates() {

		initUpdateTemplates(getTemplatesNumber());

	}

	public void initUpdateTemplates(int num) {
//		TasksManager.getInstance().addTask(new UpdateTemplatesTask(this, num));
		new UpdateTemplatesTask(this, num).execute();
		saveConfig();
	}

	public synchronized void updateTemplates(int num) {
		System.out.println("updateTemplates");
		try {
			if (getMode() != Mode.Manual)
				return;

			initTemplates(num);
			ImageGeneric igp;
			List<ManualParticle> particles;
			int size;
			ManualParticle particle;
			double[] align;

			System.err.println("JM_DEBUG: ============= Updating TEMPLATES ============");
			//FIXME: the template update needs to be done in a 
			// more efficient way, now we are using this maxcount
			// to limit the number of particles used in the update
			int count = 0;
			int maxcount = 100; 
			//FIXME: This is a tweak to avoid losing time with big particles
			if (getSize() > 256)
				maxcount = 10;
			
			for (SingleParticlePickerMicrograph m : getMicrographs()) {
				if(count >= maxcount)
					break;
				particles = m.getManualParticles();
				size = particles.size();
				//System.err.println("JM_DEBUG: Updating for Micrograph: "+ m.getFile());
				for (int i = 0; i < size; i++) {
					particle = particles.get(i);
					igp = particle.getImageGeneric();
					if (getTemplateIndex() < getTemplatesNumber())
						setTemplate(igp);
					else {
						align = templates.alignImage(igp);
						applyAlignment(particle, igp, align);
					}
					count ++;
					
				}
				saveTemplates(); //Save templates after each micrograph?
			}
			templates.getRadialAvg(radialtemplates);
			saveTemplates();
			// System.out.println("templates updated");
		} catch (Exception e) {
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public synchronized void addParticleToTemplates(ManualParticle particle) {
		try {
			ImageGeneric igp = particle.getImageGeneric();
			// will happen only in manual mode
			if (getTemplateIndex() < getTemplatesNumber())
				setTemplate(igp);
			else {
				double[] align = getTemplates().alignImage(igp);
				applyAlignment(particle, igp, align);
			}
			templates.getRadialAvg(radialtemplates);
			saveTemplates();

		} catch (Exception e) {
			throw new IllegalArgumentException(e);
		}

	}

	public boolean hasManualParticles() {
		for (SingleParticlePickerMicrograph m : micrographs) {
			if (m.hasManualParticles())
				return true;
		}
		return false;
	}

	public int getTemplateIndex() {
		return templateindex;
	}

	public synchronized void centerParticle(ManualParticle p) {
		
		if (getManualParticlesNumber() <= getTemplatesNumber())
			return;// templates not ready
		Particle shift = null;
		try {
			ImageGeneric igp = p.getImageGeneric();

			shift = radialtemplates.bestShift(igp);
			double distance = Math.sqrt(Math.pow(shift.getX(), 2)
					+ Math.pow(shift.getY(), 2))
					/ getSize();
			// System.out.printf("normalized distance:%.2f\n", distance);
			if (distance < 0.25) {
				p.setX(p.getX() + shift.getX());
				p.setY(p.getY() + shift.getY());
			}

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public synchronized void applyAlignment(ManualParticle particle,
			ImageGeneric igp, double[] align) {
		try {
			particle.setLastalign(align);
			templates.applyAlignment(igp, particle.getTemplateIndex(),
					particle.getTemplateRotation(), particle.getTemplateTilt(),
					particle.getTemplatePsi());
			// System.out.printf("adding particle: %d %.2f %.2f %.2f\n",
			// particle.getTemplateIndex(), particle.getTemplateRotation(),
			// particle.getTemplateTilt(), particle.getTemplatePsi());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void loadMicrographData(SingleParticlePickerMicrograph micrograph) {
		try {
			String file = getOutputPath(micrograph.getPosFile());
			MicrographState state;
			Integer autopickpercent;
			if (!new File(file).exists())
				return;
			if (MetaData.containsBlock(file, "header")) {
				MetaData md = new MetaData("header@" + file);
				boolean hasautopercent = md.containsLabel(MDLabel.MDL_PICKING_AUTOPICKPERCENT);
				long id = md.firstObject();
				state = MicrographState.valueOf(md.getValueString(
						MDLabel.MDL_PICKING_MICROGRAPH_STATE, id));
				if (hasautopercent)
					autopickpercent = md.getValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, id);
				else
					autopickpercent = 50;// compatibility with previous projects
				double threshold = md.getValueDouble(MDLabel.MDL_COST, id);
				micrograph.setThreshold(threshold);
				micrograph.setState(state);
				micrograph.setAutopickpercent(autopickpercent);
				md.destroy();
			}
			loadManualParticles(micrograph, file);
			loadAutomaticParticles(micrograph, file);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void loadManualParticles(SingleParticlePickerMicrograph micrograph,
			String file) {
		if (new File(file).exists()
				&& MetaData.containsBlock(file,
						getParticlesBlockName(Format.Xmipp301))) {

			int x, y;
			ManualParticle particle;

			try {
				MetaData md = new MetaData(getParticlesBlock(file));

				for (long id : md.findObjects()) {
					x = md.getValueInt(MDLabel.MDL_XCOOR, id);
					y = md.getValueInt(MDLabel.MDL_YCOOR, id);
					particle = new ManualParticle(x, y, this, micrograph);
					micrograph.addManualParticle(particle, this, false);
				}
				md.destroy();
			} catch (Exception e) {
				getLogger().log(Level.SEVERE, e.getMessage(), e);
				throw new IllegalArgumentException(e.getMessage());
			}
		}
	}

	/**
	 * This function will be called when the user selects Autopick for the first
	 * time. In order to train the classifier a metadata will be create with a
	 * list of micrographs and the corresponding .pos file associated with the
	 * micrograph, but only for those micrographs containing manual particles.
	 * 
	 * @param frame
	 * @param autopickout
	 */
	public void trainAndAutopick(SingleParticlePickerJFrame frame,
			Rectangle autopickout) {
		frame.getCanvas().setEnabled(false);
		XmippWindowUtil.blockGUI(frame, "Training and Autopicking...");
		MetaData trainmd = new MetaData();
		// Change the state of the micrograph and save changes
		micrograph.setState(MicrographState.Supervised);
		saveData(micrograph);
		// Create the metadata with each micrograph and its .pos file
		// Add all micrographs with manual particles except the current
		// micrograph
		for (SingleParticlePickerMicrograph m : micrographs) {
			if (m.hasManualParticles() && !m.equals(micrograph))
				addMicrographPos(trainmd, m);
		}
		// Add the current micrograph as the last one
		addMicrographPos(trainmd, micrograph);

		new Thread(new TrainRunnable(frame, trainmd, autopickout))
				.start();
		// new TrainRunnable(frame, trainmd, outputmd).run();
	}

	/** Helper function to add a micrograph and its .pos file */
	private void addMicrographPos(MetaData md,
			SingleParticlePickerMicrograph micrograph) {
		String posfile = getOutputPath(micrograph.getPosFile());
		long id = md.addObject();
		md.setValueString(MDLabel.MDL_MICROGRAPH, micrograph.getFile(), id);
		md.setValueString(MDLabel.MDL_MICROGRAPH_PARTICLES, posfile, id);

	}

	public class TrainRunnable implements Runnable {

		private SingleParticlePickerJFrame frame;
		private MetaData trainmd;
		private MetaData outputmd;
		private Rectangle autopickout;

		public TrainRunnable(SingleParticlePickerJFrame frame,
				MetaData trainmd, Rectangle autopickout) {
			this.frame = frame;
			this.trainmd = trainmd;
			this.outputmd = new MetaData();
			this.autopickout = autopickout;
		}

		public void run() {
			try {
				classifier.train(trainmd, (int) autopickout.getX(),
						(int) autopickout.getY(), (int) autopickout.getWidth(),
						(int) autopickout.getHeight());// should remove training
														// files
				micrograph.setAutopickpercent(autopickpercent);
				classifier.autopick(micrograph.getFile(), outputmd,
						micrograph.getAutopickpercent());
				int x, y;
				double cost;
				AutomaticParticle ap;
				for (long id : outputmd.findObjects()) {
					x = outputmd.getValueInt(MDLabel.MDL_XCOOR, id);
					y = outputmd.getValueInt(MDLabel.MDL_YCOOR, id);
					cost = outputmd.getValueDouble(MDLabel.MDL_COST, id);
					ap = new AutomaticParticle(x, y, SingleParticlePicker.this,
							micrograph, cost, false);
					if (autopickout.contains(new Point(x, y))) {
						ap.setDeleted(true);
						ap.setCost(-1);
					}
					micrograph.addAutomaticParticle(ap);
				}

				XmippWindowUtil.releaseGUI(frame.getRootPane());
				frame.getCanvas().setEnabled(true);
				frame.getCanvas().repaint();
				frame.updateMicrographsModel();
				
				trainmd.destroy();
				outputmd.destroy();
			} catch (Exception e) {
				ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
				throw new IllegalArgumentException(e.getMessage());
			}
		}
	}

	public void autopick(SingleParticlePickerJFrame frame,
			SingleParticlePickerMicrograph next) {
		next.setState(MicrographState.Supervised);
		saveData(next);

		frame.getCanvas().setEnabled(false);
		XmippWindowUtil.blockGUI(frame, "Autopicking...");

		new Thread(new AutopickRunnable(frame, next)).start();

	}

	public class AutopickRunnable implements Runnable {

		private SingleParticlePickerJFrame frame;
		private MetaData outputmd;
		private SingleParticlePickerMicrograph micrograph;

		public AutopickRunnable(SingleParticlePickerJFrame frame,
				SingleParticlePickerMicrograph micrograph) {
			this.frame = frame;
			this.outputmd = new MetaData();
			this.micrograph = micrograph;
		}

		public void run() {
			micrograph.getAutomaticParticles().clear();
			micrograph.setAutopickpercent(autopickpercent);
			classifier.autopick(micrograph.getFile(), outputmd,
					micrograph.getAutopickpercent());
			loadAutomaticParticles(micrograph, outputmd);
			String path = getParticlesAutoBlock(micrograph);
			outputmd.writeBlock(path);
			frame.getCanvas().repaint();
			frame.getCanvas().setEnabled(true);
			frame.updateMicrographsModel();
			XmippWindowUtil.releaseGUI(frame.getRootPane());
			outputmd.destroy();
		}

	}

	public void correctAndAutopick(SingleParticlePickerJFrame frame,
			SingleParticlePickerMicrograph current,
			SingleParticlePickerMicrograph next, Rectangle correctout) {
		current.setState(MicrographState.Corrected);
		if (getMode() == Mode.Supervised
				&& next.getState() == MicrographState.Available)
			next.setState(MicrographState.Supervised);
		saveData(current);
		saveData(next);
		MetaData addedmd = getAddedMetaData(current, correctout);
		MetaData automd = new MetaData(getParticlesAutoBlock(micrograph));
		frame.getCanvas().setEnabled(false);
		XmippWindowUtil.blockGUI(frame, "Correcting and Autopicking...");
		new Thread(new CorrectAndAutopickRunnable(frame, addedmd, automd, next)).start();

	}

	private MetaData getAddedMetaData(SingleParticlePickerMicrograph m,
			Rectangle correctout) {
		MetaData addedmd = null;
		if (correctout == null)
			addedmd = new MetaData(getParticlesBlock(micrograph));
		else {
			long id;
			addedmd = new MetaData();
			for (ManualParticle p : m.getManualParticles())
				if (!correctout.contains(new Point(p.getX(), p.getY()))) {
					id = addedmd.addObject();
					addedmd.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
					addedmd.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
				}
		}
		return addedmd;
	}

	public class CorrectAndAutopickRunnable implements Runnable {

		private MetaData manualmd;
		private MetaData automaticmd;
		private SingleParticlePickerMicrograph next;
		private SingleParticlePickerJFrame frame;
		private MetaData outputmd;

		public CorrectAndAutopickRunnable(SingleParticlePickerJFrame frame,
				MetaData manualmd, MetaData automaticmd,
				SingleParticlePickerMicrograph next) {
			this.frame = frame;
			this.manualmd = manualmd;
			this.automaticmd = automaticmd;
			this.next = next;
			this.outputmd = new MetaData();
		}

		public void run() {

			classifier.correct(manualmd, automaticmd, micrograph.getThreshold());
			manualmd.destroy();
			automaticmd.destroy();
			
			if (getMode() == Mode.Supervised
					&& next.getState() == MicrographState.Supervised) {
				next.getAutomaticParticles().clear();
				next.setAutopickpercent(autopickpercent);
				classifier.autopick(next.getFile(), outputmd,
						next.getAutopickpercent());
				loadAutomaticParticles(next, outputmd);
				String path = getParticlesAutoBlock(next);
				outputmd.writeBlock(path);
			}
			outputmd.destroy();
			frame.getCanvas().repaint();
			frame.getCanvas().setEnabled(true);
			frame.updateMicrographsModel(true);
			XmippWindowUtil.releaseGUI(frame.getRootPane());

		}
	}

	public void resetMicrograph(SingleParticlePickerMicrograph micrograph) {
		micrograph.getManualParticles().clear();
		micrograph.getAutomaticParticles().clear();
		micrograph.setState(MicrographState.Available);
		micrograph.setThreshold(0);
		new File(getOutputPath(micrograph.getPosFile())).delete();

	}

	public void correct(Rectangle correctout) {
		micrograph.setState(MicrographState.Corrected);
		saveData(micrograph);
		MetaData addedmd = getAddedMetaData(micrograph, correctout);
		MetaData automd = new MetaData(getParticlesAutoBlock(micrograph));
		classifier.correct(addedmd, automd, micrograph.getThreshold());
		addedmd.destroy();
		automd.destroy();
	}

}
