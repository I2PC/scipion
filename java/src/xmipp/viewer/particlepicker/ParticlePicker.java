package xmipp.viewer.particlepicker;

import ij.CommandListener;
import ij.Executer;
import ij.ImageListener;
import ij.ImagePlus;
import ij.plugin.frame.Recorder;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.jni.Particle;
import xmipp.jni.Program;
import xmipp.utils.TasksManager;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;
import xmipp.viewer.particlepicker.training.model.TrainingPicker;

public abstract class ParticlePicker {

	protected String familiesfile;
	protected String macrosfile;
	protected static Logger logger;
	protected String outputdir = ".";
	protected boolean changed;
	protected List<Family> families;
	protected FamilyState mode;
	protected List<IJCommand> filters;
	protected String selfile;
	protected String command;
	protected Family family;
	protected String configfile;
	public static final int defAutoPickPercent = 90;
	protected int autopickpercent = defAutoPickPercent;
	
	public static final int fsizemax = 800;
	protected String block;
	
	public ParticlePicker(String selfile, String outputdir, FamilyState mode) {
		this(selfile, outputdir, null, mode);

	}

	public ParticlePicker(String selfile, String outputdir, String fname, FamilyState mode)
	{
		this(null, selfile, outputdir, fname, mode);
	}

	public ParticlePicker(String block, String selfile, String outputdir, String fname, FamilyState mode)
	{
		this.block = block;
		this.outputdir = outputdir;
		this.selfile = selfile;
		this.mode = mode;
		this.configfile = getOutputPath("config.xmd");
		initFamilies(fname);
		initFilters();
		loadEmptyMicrographs();
		loadConfig();
	}
	
	
	
	protected void initFamilies(String fname)
	{
		this.familiesfile = getOutputPath("families.xmd");
		
		this.families = new ArrayList<Family>();
		
		loadFamilies();
		if (fname == null)
			family = families.get(0);
		else
			family = getFamily(fname);
		if (family == null)
			throw new IllegalArgumentException("Invalid family " + fname);

	}
	
	public int getSize() {
		return family.getSize();
	}

	public Color getColor() {
		return family.getColor();
	}

	public void setColor(Color color) {
		family.setColor(color);
		saveFamilies();
	}

	public void setSize(int size) {
		family.setSize(size);
		saveFamilies();
	}

	public Family getFamily() {
		return family;
	}

	public void setFamily(Family family) {
		this.family = family;

	}

	public String getPosFileFromXmipp24Project(String projectdir, String mname) {
		String suffix = ".raw.Common.pos";
		return String.format("%1$s%2$sPreprocessing%2$s%3$s%2$s%3$s%4$s", projectdir, File.separator, mname, suffix);
	}

	
	public abstract void loadEmptyMicrographs();

	private void initFilters() {
		this.macrosfile = getOutputPath("macros.xmd");
		filters = new ArrayList<IJCommand>();
		loadFilters();
		Recorder.record = true;

		// detecting if a command is thrown by ImageJ
		Executer.addCommandListener(new CommandListener() {
			public String commandExecuting(String command) {
				ParticlePicker.this.command = command;
				return command;

			}
		});
		ImagePlus.addImageListener(new ImageListener() {

			@Override
			public void imageUpdated(ImagePlus arg0) {
				updateFilters();
			}

			@Override
			public void imageOpened(ImagePlus arg0) {

			}

			@Override
			public void imageClosed(ImagePlus arg0) {
				// TODO Auto-generated method stub

			}
		});
	}

	private void updateFilters() {
		if (command != null) {
			String options = "";
			if (Recorder.getCommandOptions() != null) options = Recorder.getCommandOptions();
			if (!isFilterSelected(command))
				addFilter(command, options);
			else if (!(options == null || options.equals(""))) for (IJCommand f : filters)
				if (f.getCommand().equals(command)) f.setOptions(options);
			saveFilters();
			command = null;

		}
	}

	public String getMicrographsSelFile() {
		return selfile;
	}

	public void addFilter(String command, String options) {
		IJCommand f = new IJCommand(command, options);
		filters.add(f);
	}

	public List<IJCommand> getFilters() {
		return filters;
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
//				FileHandler fh = new FileHandler("PPicker.log", true);
//				fh.setFormatter(new SimpleFormatter());
				logger = Logger.getLogger("PPickerLogger");
//				logger.addHandler(fh);
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

	public String getOutputDir() {
		return outputdir;
	}

	public List<Family> getFamilies() {
		return families;
	}

	public void saveFamilies() {
		long id;
		String file = familiesfile;
		try {
			MetaData md = new MetaData();
			for (Family f : families) {
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, f.getName(), id);
				md.setValueInt(MDLabel.MDL_COLOR, f.getColor().getRGB(), id);
				md.setValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, f.getSize(), id);
				md.setValueInt(MDLabel.MDL_PICKING_FAMILY_TEMPLATES, f.getTemplatesNumber(), id);//saved for compatibility with future versions
				md.setValueString(MDLabel.MDL_PICKING_FAMILY_STATE, f.getStep().toString(), id);
			}
			md.write(file);
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public abstract List<? extends Micrograph> getMicrographs();
	
	public int getMicrographIndex()
	{
		return getMicrographs().indexOf(getMicrograph());
	}

	public void loadFamilies() {
		families.clear();
		String file = familiesfile;
		if (!new File(file).exists()) {
			families.add(new Family("DefaultFamily", Color.green, fsizemax/4, 1, getTemplatesFile("DefaultFamily")));
			saveFamilies();
			return;
		}

		Family family;
		int rgb, size;
		Integer templatesNumber = 1;
		FamilyState state;
		String name, templatesfile;
		ImageGeneric templates;
		try {
			MetaData md = new MetaData(file);
			long[] ids = md.findObjects();
			for (long id : ids) {
				name = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				rgb = md.getValueInt(MDLabel.MDL_COLOR, id);
				size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
				templatesNumber = md.getValueInt(MDLabel.MDL_PICKING_FAMILY_TEMPLATES, id);
				if( templatesNumber == null || templatesNumber == 0)
					templatesNumber = 1;//for compatibility with previous projects
				state = FamilyState.valueOf(md.getValueString(MDLabel.MDL_PICKING_FAMILY_STATE, id));
				state = validateState(state);
				templatesfile = getTemplatesFile(name);
				if (new File(templatesfile).exists() )
				{
					templates = new ImageGeneric(templatesfile);
					family = new Family(name, new Color(rgb), size, state, this, templates);
					
				}
				else
				{
					family = new Family(name, new Color(rgb), size, state, this, templatesNumber, templatesfile);
				}
				families.add(family);
			}
			md.destroy();
			if (families.size() == 0) throw new IllegalArgumentException(String.format("No families specified on %s", file));
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public String getTemplatesFile(String name)
	{
		return getOutputPath(name + "_templates.stk");
	}

	public FamilyState validateState(FamilyState state) {

		if (mode == FamilyState.Review && state != FamilyState.Review) {
			setChanged(true);
			return FamilyState.Review;
		}
		if (mode == FamilyState.Manual && !(state == FamilyState.Manual || state == FamilyState.Available))
			throw new IllegalArgumentException(String.format("Can not use %s mode on this data", mode));
		if (mode == FamilyState.Supervised && state == FamilyState.Review)
			throw new IllegalArgumentException(String.format("Can not use %s mode on this data", mode));
		return state;

	}// function validateState

	public Family getFamily(String name) {
		if(name == null)
			return null;
		for (Family f : getFamilies())
			if (f.getName().equalsIgnoreCase(name)) return f;
		return null;
	}// function getFamily

	public boolean existsFamilyName(String name) {
		return getFamily(name) != null;
	}// function existsFamilyName

	protected boolean containsBlock(String file, String block) {
		try {
			return Arrays.asList(MetaData.getBlocksInMetaDataFile(file)).contains(block);
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}// function containsBlock

	public void saveData() {
		saveFilters();
		saveFamilies();
	}// function saveData

	public abstract void saveData(Micrograph m);

	public abstract int getManualParticlesNumber(Family f);

	public void saveFilters() {
		long id;
		String file = macrosfile;
		if (filters.isEmpty()) {
			new File(file).delete();
			return;
		}
		String options;
		try {
			MetaData md = new MetaData();
			for (IJCommand f : filters) {
				id = md.addObject();
				md.setValueString(MDLabel.MDL_IMAGE1, f.getCommand().replace(' ', '_'), id);
				options = (f.getOptions() == null || f.getOptions().equals("")) ? "NULL" : f.getOptions().replace(' ', '_');
				md.setValueString(MDLabel.MDL_IMAGE2, options, id);
			}
			md.write(file);
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}// function persistFilters

	public void loadFilters() {
		filters.clear();
		String file = macrosfile;
		if (!new File(file).exists()) return;

		String command, options;
		try {
			MetaData md = new MetaData(file);
			long[] ids = md.findObjects();
			for (long id : ids) {
				command = md.getValueString(MDLabel.MDL_IMAGE1, id).replace('_', ' ');
				options = md.getValueString(MDLabel.MDL_IMAGE2, id).replace('_', ' ');
				if (options.equals("NULL")) options = "";
				filters.add(new IJCommand(command, options));

			}
			md.destroy();
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}// function loadFilters

	public Micrograph getMicrograph(String name) {
		for (Micrograph m : getMicrographs())
			if (m.getName().equalsIgnoreCase(name)) return m;
		return null;
	}

	public void setAutopickpercent(int autopickpercent)
	{
		this.autopickpercent = autopickpercent;
	}
	
	public int getAutopickpercent()
	{
		return autopickpercent;
	}
	
	

	public abstract void loadConfig();
	
	public abstract void saveConfig();
	

	void removeFilter(String filter) {
		for (IJCommand f : filters)
			if (f.getCommand().equals(filter)) {
				filters.remove(f);
				saveFilters();
				break;
			}
	}// function removeFilter

	public boolean isFilterSelected(String filter) {
		for (IJCommand f : filters)
			if (f.getCommand().equals(filter)) return true;
		return false;
	}// function isFilterSelected

	public abstract void exportParticles(String absolutePath);

	public Format detectFormat(String path) {
		return Format.Unknown;
	}

	public Format detectFileFormat(String path) {
		if (path.endsWith(".raw.Common.pos")) return Format.Xmipp24;
		if (path.endsWith(".pos")) return Format.Xmipp30;
		if (path.endsWith(".box")) return Format.Eman;
		return Format.Unknown;
	}

	public abstract String getImportMicrographName(String path, String filename, Format f);

	/** Return the number of particles imported */
	public abstract String importParticlesFromFolder(String path, Format f, float scale, boolean invertx, boolean inverty);

	/** Return the number of particles imported from a file */
	public void fillParticlesMdFromFile(String path, Format f, Micrograph m, MetaData md, float scale, boolean invertx, boolean inverty) {

		if (f == Format.Auto) f = detectFileFormat(path);

		switch (f) {
		case Xmipp24:
			md.readPlain(path, "xcoor ycoor");
			break;
		case Xmipp30:
			if(!containsBlock(path, family.getName()))
				throw new IllegalArgumentException(XmippMessage.getIllegalValueMsgWithInfo("family", family.getName(), "Particles for this family are not defined in file"));
			md.read(String.format("%s@%s", family.getName(), path));
			break;
		case Eman:
			fillParticlesMdFromEmanFile(path, m, md, scale);
			break;
		default:
			md.clear();
		}
		int width = (int) (m.width / scale);// original width
		int height = (int) (m.height / scale);// original height
		if (invertx) md.operate(String.format("xcoor=%d-xcoor", width));
		if (inverty) md.operate(String.format("ycoor=%d-ycoor", height));
		if (scale != 1.f) md.operate(String.format("xcoor=xcoor*%f,ycoor=ycoor*%f", scale, scale));

	}// function importParticlesFromFile
	
	

	public void fillParticlesMdFromEmanFile(String file, Micrograph m, MetaData md, float scale) {
		String line = "";
		// System.out.println("Importing from EMAN, file: " + file);
		try {
			BufferedReader reader = null;
			reader = new BufferedReader(new FileReader(file));
			line = reader.readLine();
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		// inverty = true;
		md.readPlain(file, "xcoor ycoor particleSize");

		long fid = md.firstObject();
		int size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, fid);
		if (size > 0) family.setSize(Math.round(size * scale));
		int half = size / 2;
		md.operate(String.format("xcoor=xcoor+%d,ycoor=ycoor+%d", half, half));

	}// function fillParticlesMdFromEmanFile

	public void runXmippProgram(String program, String args) {
		try {
			Program.runByName(program, args);
		} catch (Exception e) {
			TrainingPicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	public abstract Micrograph getMicrograph();

	public abstract void setMicrograph(Micrograph m);
	
		


	
	
	

	public abstract boolean isValidSize(int size);


}
