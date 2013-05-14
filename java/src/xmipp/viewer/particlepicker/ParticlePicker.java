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
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.jni.Program;
import xmipp.viewer.particlepicker.training.model.Mode;

public abstract class ParticlePicker
{

	protected String macrosfile;
	protected static Logger logger;
	protected String outputdir = ".";
	protected boolean changed;
	protected Mode mode;
	protected List<IJCommand> filters;
	protected String selfile;
	protected String command;
	protected String configfile;
	public static final int defAutoPickPercent = 90;
	protected int autopickpercent = defAutoPickPercent;
	
	private Color color;
	private int size;

	
	

	public static final int sizemax = 800;
	protected String block;
	
	
	private static Color[] colors = new Color[] { Color.BLUE, Color.CYAN,
		Color.GREEN, Color.MAGENTA, Color.ORANGE, Color.PINK, Color.YELLOW };
	
	private static int nextcolor;
	
	public static Color getNextColor() {
		Color next = colors[nextcolor];
		nextcolor++;
		if (nextcolor == colors.length)
			nextcolor = 0;
		return next;
	}
	
	
	
	
	public ParticlePicker(String selfile, Mode mode)
	{
		this(null, selfile, ".", mode);
	}
	
	public ParticlePicker(String selfile, String outputdir, Mode mode)
	{
		this(null, selfile, outputdir, mode);
	}



	public ParticlePicker(String block, String selfile, String outputdir, Mode mode)
	{
		this.block = block;
		this.outputdir = outputdir;
		this.selfile = selfile;
		this.outputdir = outputdir;
		this.mode = mode;
		this.configfile = getOutputPath("config.xmd");
		color = getNextColor();
		size = getDefaultSize();
		initFilters();
		loadEmptyMicrographs();
		loadConfig();
	}
	
	

	public void loadConfig()
	{
		String file = configfile;
		if (!new File(file).exists())
		{
			setMicrograph(getMicrographs().get(0));
			return;

		}

		String mname;
		int rgb;
		try
		{
			MetaData md = new MetaData(file);
			for (long id : md.findObjects())
			{
		
				mname = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
				setMicrograph(getMicrograph(mname));
				rgb = md.getValueInt(MDLabel.MDL_COLOR, id);
				color = new Color(rgb);
				size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}





	public String getPosFileFromXmipp24Project(String projectdir, String mname)
	{
		String suffix = ".raw.Common.pos";
		return String.format("%1$s%2$sPreprocessing%2$s%3$s%2$s%3$s%4$s", projectdir, File.separator, mname, suffix);
	}

	public abstract void loadEmptyMicrographs();

	private void initFilters()
	{

		this.macrosfile = getOutputPath("macros.xmd");
		filters = new ArrayList<IJCommand>();
		loadFilters();
		Recorder.record = true;

		// detecting if a command is thrown by ImageJ
		Executer.addCommandListener(new CommandListener()
		{
			public String commandExecuting(String command)
			{
				ParticlePicker.this.command = command;
				return command;

			}
		});
		ImagePlus.addImageListener(new ImageListener()
		{

			@Override
			public void imageUpdated(ImagePlus arg0)
			{
				updateFilters();
			}

			@Override
			public void imageOpened(ImagePlus arg0)
			{

			}

			@Override
			public void imageClosed(ImagePlus arg0)
			{
				// TODO Auto-generated method stub

			}
		});
	}

	private void updateFilters()
	{
		if (command != null)
		{
			String options = "";
			if (Recorder.getCommandOptions() != null)
				options = Recorder.getCommandOptions();
			if (!isFilterSelected(command))
				addFilter(command, options);
			else if (!(options == null || options.equals("")))
				for (IJCommand f : filters)
					if (f.getCommand().equals(command))
						f.setOptions(options);

			saveFilters();
			command = null;

		}
	}

	public String getMicrographsSelFile()
	{
		return selfile;
	}

	public void addFilter(String command, String options)
	{
		IJCommand f = new IJCommand(command, options);
		filters.add(f);
	}

	public List<IJCommand> getFilters()
	{
		return filters;
	}

	public void setChanged(boolean changed)
	{
		this.changed = changed;
	}

	public boolean isChanged()
	{
		return changed;
	}

	public Mode getMode()
	{
		return mode;
	}

	public static Logger getLogger()
	{
		try
		{
			if (logger == null)
			{
				FileHandler fh = new FileHandler("PPicker.log", true);
				fh.setFormatter(new SimpleFormatter());
				logger = Logger.getLogger("PPickerLogger");
				// logger.addHandler(fh);
			}
			return logger;
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	public String getOutputPath(String file)
	{
		return outputdir + File.separator + file;
	}

	public String getOutputDir()
	{
		return outputdir;
	}



	public abstract List<? extends Micrograph> getMicrographs();

	public int getMicrographIndex()
	{
		return getMicrographs().indexOf(getMicrograph());
	}



	public String getTemplatesFile(String name)
	{
		return getOutputPath(name + "_templates.stk");
	}

	public Mode validateState(Mode state) {


		if (mode == Mode.Review && state != Mode.Review)
		{
			setChanged(true);
			return Mode.Review;
		}
		if (mode == Mode.Manual && !(state == Mode.Manual || state == Mode.Available))
			throw new IllegalArgumentException(String.format("Can not use %s mode on this data", mode));
		if (mode == Mode.Supervised && state == Mode.Review)
			throw new IllegalArgumentException(String.format("Can not use %s mode on this data", mode));
		return state;

	}// function validateState


	public void saveData()
	{
		saveFilters();
		saveConfig();
	}// function saveData

	public abstract void saveData(Micrograph m);

	


	public void saveFilters() {
		long id;
		String file = macrosfile;
		if (filters.isEmpty())
		{
			new File(file).delete();
			return;
		}
		String options;
		try
		{
			MetaData md = new MetaData();
			for (IJCommand f : filters)
			{
				id = md.addObject();
				md.setValueString(MDLabel.MDL_IMAGE1, f.getCommand().replace(' ', '_'), id);
				options = (f.getOptions() == null || f.getOptions().equals("")) ? "NULL" : f.getOptions().replace(' ', '_');
				md.setValueString(MDLabel.MDL_IMAGE2, options, id);
			}
			md.write(file);
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}// function persistFilters

	public void loadFilters()
	{
		filters.clear();
		String file = macrosfile;
		if (!new File(file).exists())
			return;

		String command, options;
		try
		{
			MetaData md = new MetaData(file);
			long[] ids = md.findObjects();
			for (long id : ids)
			{
				command = md.getValueString(MDLabel.MDL_IMAGE1, id).replace('_', ' ');
				options = md.getValueString(MDLabel.MDL_IMAGE2, id).replace('_', ' ');
				if (options.equals("NULL"))
					options = "";
				filters.add(new IJCommand(command, options));

			}
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}// function loadFilters

	public Micrograph getMicrograph(String name)
	{
		for (Micrograph m : getMicrographs())
			if (m.getName().equalsIgnoreCase(name))
				return m;
		return null;
	}



	public void saveConfig()
	{
		try
		{
			MetaData md;
			String file = configfile;
			md = new MetaData();
			long id = md.addObject();
			saveConfig(md, id);
			md.write(file);
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
		md.setValueString(MDLabel.MDL_MICROGRAPH, getMicrograph().getName(), id);
		md.setValueInt(MDLabel.MDL_COLOR, getColor().getRGB(), id);
		md.setValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, getSize(), id);
		
	}

	void removeFilter(String filter)
	{
		for (IJCommand f : filters)
			if (f.getCommand().equals(filter))
			{
				filters.remove(f);
				saveFilters();
				break;
			}
	}// function removeFilter

	public boolean isFilterSelected(String filter)
	{
		for (IJCommand f : filters)
			if (f.getCommand().equals(filter))
				return true;
		return false;
	}// function isFilterSelected

	

	public Format detectFormat(String path)
	{
		return Format.Unknown;
	}

	public Format detectFileFormat(String path)
	{
		if (path.endsWith(".raw.Common.pos"))
			return Format.Xmipp24;
		if (path.endsWith(".pos"))
			return Format.Xmipp30;
		if (path.endsWith(".box"))
			return Format.Eman;
		return Format.Unknown;
	}

	/** Return the number of particles imported from a file */
	public void fillParticlesMdFromFile(String path, Format f, Micrograph m, MetaData md, float scale, boolean invertx, boolean inverty)
	{

		if (f == Format.Auto)
			f = detectFileFormat(path);

		switch (f)
		{
		case Xmipp24:
			md.readPlain(path, "xcoor ycoor");
			break;
//		case Xmipp30:
//			if (!MetaData.containsBlock(path, family.getName()))
//				throw new IllegalArgumentException(
//						XmippMessage.getIllegalValueMsgWithInfo("family", family.getName(), "Particles for this family are not defined in file"));
//			md.read(String.format("%s@%s", family.getName(), path));
//			break;
		case Eman:
			fillParticlesMdFromEmanFile(path, m, md, scale);
			break;
		default:
			md.clear();
		}
		int width = (int) (m.width / scale);// original width
		int height = (int) (m.height / scale);// original height
		if (invertx)
			md.operate(String.format("xcoor=%d-xcoor", width));
		if (inverty)
			md.operate(String.format("ycoor=%d-ycoor", height));
		if (scale != 1.f)
			md.operate(String.format("xcoor=xcoor*%f,ycoor=ycoor*%f", scale, scale));

	}// function importParticlesFromFile

	public void fillParticlesMdFromEmanFile(String file, Micrograph m, MetaData md, float scale)
	{
		String line = "";
		// System.out.println("Importing from EMAN, file: " + file);
		try
		{
			BufferedReader reader = null;
			reader = new BufferedReader(new FileReader(file));
			line = reader.readLine();
			reader.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		// inverty = true;
		md.readPlain(file, "xcoor ycoor particleSize");

		long fid = md.firstObject();
		int size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, fid);
		if (size > 0)
//			family.setSize(Math.round(size * scale));
			setSize(Math.round(size * scale));
		int half = size / 2;
		md.operate(String.format("xcoor=xcoor+%d,ycoor=ycoor+%d", half, half));

	}// function fillParticlesMdFromEmanFile


	

	public void runXmippProgram(String program, String args)
	{
		try
		{
			Program.runByName(program, args);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	public abstract Micrograph getMicrograph();

	public abstract void setMicrograph(Micrograph m);
	
		

	
	public abstract boolean isValidSize(int size);
	
	
	
	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		if (size >  ParticlePicker.sizemax)
			throw new IllegalArgumentException(String.format(
					"Max size is %s, %s not allowed",  ParticlePicker.sizemax, size));

		this.size = size;
		saveConfig();//size and color will be on config
	}
	


	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
		saveConfig();
	}
	
	public static Color getColor(String name) {
		Color color;
		try {
			Field field = Class.forName("java.awt.Color").getField(name);
			color = (Color) field.get(null);// static field, null for parameter
		} catch (Exception e) {
			color = null; // Not defined
		}
		return color;
	}

	public static int getDefaultSize() {
		return 100;
	}



	public int getRadius() {
		return size / 2;
	}


}
