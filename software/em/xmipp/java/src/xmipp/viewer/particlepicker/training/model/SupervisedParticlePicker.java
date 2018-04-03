package xmipp.viewer.particlepicker.training.model;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;

import javax.swing.JFrame;
import javax.swing.SwingWorker;

import xmipp.jni.Classifier;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MDRow;
import xmipp.jni.MetaData;
import xmipp.jni.Particle;
import xmipp.jni.PickingClassifier;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.JMetaDataIO;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.particlepicker.*;
import xmipp.viewer.particlepicker.training.AutopickRunnable;
import xmipp.viewer.particlepicker.training.CorrectAndAutopickRunnable;
import xmipp.viewer.particlepicker.training.TrainRunnable;
import xmipp.viewer.particlepicker.training.gui.SupervisedPickerJFrame;
import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;
import xmipp.viewer.scipion.ScipionMetaData;

/**
 * Business object for Single Particle Picker GUI. Inherits from ParticlePicker
 * 
 * @author airen
 * 
 */
public class SupervisedParticlePicker extends ParticlePicker
{

	protected List<SupervisedPickerMicrograph> micrographs;
	private SupervisedPickerMicrograph micrograph;
	protected int autopickpercent;

        //used in previous versions
	private Integer threads = 1;
	private boolean fastmode = true;
	private boolean incore = false;
        ///////////////////////////

	public static int dtemplatesnum = 1;
	private ImageGeneric templates;
	private ImageGeneric radialtemplates;
	private String templatesfile;
	private int templateindex;
	private Classifier classifier;
    private UpdateTemplatesTask uttask;
    private TemplatesJDialog dialog;
    private boolean isautopick;
    private boolean existspsd;
    private boolean hasDefocusU;

	// private String reviewfile;

	public SupervisedParticlePicker(String selfile, String outputdir, ParticlePickerParams params)
	{
		this(null, selfile, outputdir, params);

	}

	public SupervisedParticlePicker(String block, String selfile, String outputdir, ParticlePickerParams params)
	{

		super(block, selfile, outputdir, params);
		try
		{
			templatesfile = getOutputPath("templates.stk");
			if (!new File(templatesfile).exists())
				initTemplates(dtemplatesnum);
			else
			{
				this.templates = new ImageGeneric(templatesfile);
				templates.read(templatesfile, false);
				radialtemplates = new ImageGeneric(ImageGeneric.Float);
				radialtemplates.resize(getSize(), getSize(), 1, getTemplatesNumber());
                templateindex = (templates.getStatistics()[2] == 0)? 0: getTemplatesNumber();
			}
            templates.getRadialAvg(radialtemplates);
            MDRow[] micsmd = new MDRow[micrographs.size()];
            MDRow row; int index = 0;

			if(params.classifierProperties != null)
			{
				classifier = new GenericClassifier(params.classifierProperties);
				for (SupervisedPickerMicrograph m : micrographs)
					loadMicrographData(m);
			}
			else
			{
				for (SupervisedPickerMicrograph m : micrographs) {
					loadMicrographData(m);
					row = new MDRow();
					row.setValueString(MDLabel.MDL_MICROGRAPH, m.getFile());
					micsmd[index] = row;
					index ++;
				}
				classifier = new PickingClassifier(getSize(), getOutputPath("model"), micsmd);
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	
	public SupervisedParticlePicker(String selfile, String outputdir, Integer threads, boolean fastmode, boolean incore, ParticlePickerParams params)
	{
		this(selfile, outputdir, params);

		this.threads = threads;
		this.fastmode = fastmode;
		this.incore = incore;

	}

	/* Save changes to the current micrograph(To be removed). */
	public void saveData()
	{
		super.saveData();
		saveData(micrograph);
		setChanged(false);
        setSaved(true);
	}

	/* Save changes to ALL micrographs, mainly used when importing. */
	public void saveAllData()
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

	public int getTemplatesNumber()
	{
		if (templates == null)
			return dtemplatesnum;
		try
		{
			return (int) templates.getNDim();
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return dtemplatesnum;
	}
        
	public void setSize(int size)
	{

		super.setSize(size);
		classifier.setSize(size);
	}
        
	public synchronized void initTemplates()
	{
		initTemplates(getTemplatesNumber());
	}

	public synchronized void initTemplates(int num)
	{
		if (num == 0)
			return;
		try
		{
			this.templates = new ImageGeneric(ImageGeneric.Float);
			templates.resize(getSize(), getSize(), 1, num);
			templates.write(templatesfile);
			templateindex = 0;
			radialtemplates = new ImageGeneric(ImageGeneric.Float);
			radialtemplates.resize(getSize(), getSize(), 1, num);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public synchronized void setTemplatesNumber(int num)
	{
		if (num <= 0)
			throw new IllegalArgumentException(
					XmippMessage.getIllegalValueMsgWithInfo("Templates Number", Integer.valueOf(num), "Must have at least one template"));
		updateTemplatesStack(num, false);
                saveConfig();
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
			// TODO getArrayFloat and setArrayFloat must be call from C both in
			// one function
			matrix = ig.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.FIRST_SLICE);
			templates.setArrayFloat(matrix, ImageGeneric.FIRST_IMAGE + templateindex, ImageGeneric.FIRST_SLICE);
			templateindex++;

		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	// to update templates with the right particles
	public synchronized void resetParticleImages()
	{
		for (SupervisedPickerMicrograph m : micrographs)
		{
			for (ManualParticle p : m.getManualParticles())
				p.resetImage();

		}
	}

	public synchronized void updateTemplates(ImageGeneric templates, int templateindex)
	{
                
		try
		{
			this.templates = templates;
                        radialtemplates = new ImageGeneric(ImageGeneric.Float);
			radialtemplates.resize(getSize(), getSize(), 1, getTemplatesNumber());
			templates.getRadialAvg(radialtemplates);
			this.templateindex = templateindex;
                        saveTemplates();
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public synchronized void addParticleToTemplates(ManualParticle particle)
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
			templates.getRadialAvg(radialtemplates);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}

	}
        
	public synchronized void centerParticle(ManualParticle p)
	{

		if (getManualParticlesNumber() <= getTemplatesNumber())
			return;// templates not ready
		Particle shift = null;
		try
		{
			ImageGeneric igp = p.getImageGeneric();

			shift = radialtemplates.bestShift(igp);
			double distance = Math.sqrt(Math.pow(shift.getX(), 2) + Math.pow(shift.getY(), 2)) / getSize();
			// System.out.printf("normalized distance:%.2f\n", distance);
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

	public synchronized void applyAlignment(ManualParticle particle, ImageGeneric igp, double[] align)
	{
		try
		{
			particle.setLastalign(align);
			templates.applyAlignment(igp, particle.getTemplateIndex(), particle.getTemplateRotation(), particle.getTemplateTilt(), particle
					.getTemplatePsi());
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
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
            if(templateindex == getTemplatesNumber())
                templates.write(getTemplatesFile());
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}
        
    public boolean containsPSD()
    {
        return existspsd;
    }

    public boolean hasDefocusU() {
	    return hasDefocusU;
    }

    @Override
	public void loadEmptyMicrographs()
	{
    	String selfile = getMicrographsSelFile();
		if(selfile.endsWith(".sqlite"))
			loadEmptyMicrographsFromSqlite();
		else
			loadEmptyMicrographsFromMd();
	}
    
    public void loadEmptyMicrographsFromMd()
    {
    	if (micrographs == null)
			micrographs = new ArrayList<SupervisedPickerMicrograph>();
		else
			micrographs.clear();
		SupervisedPickerMicrograph micrograph;
		String psd = null, ctf = null, filename;
		try
		{
			String selfile = getMicrographsSelFile();
			MetaData md = new MetaData(selfile);
                        
			md.removeDisabled();
                        
			int fileLabel;
			if (md.containsLabel(MDLabel.MDL_MICROGRAPH))
				fileLabel = MDLabel.MDL_MICROGRAPH;
			else if (md.containsLabel(MDLabel.MDL_IMAGE))
				fileLabel = MDLabel.MDL_IMAGE;
			else
				throw new IllegalArgumentException(String.format("Labels MDL_MICROGRAPH or MDL_IMAGE not found in metadata %s", selfile));
			existspsd = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			boolean existsctf = md.containsLabel(MDLabel.MDL_CTF_MODEL);
			hasDefocusU = md.containsLabel(MDLabel.MDL_CTF_DEFOCUSV);

			boolean hasctfinfo = existsctf || existspsd || hasDefocusU;

			long[] ids = md.findObjects();
			CtfInfo ctfInfo;

			for (long id : ids)
			{
				ctfInfo = hasctfinfo ? new CtfInfo() : null;

				filename = md.getValueString(fileLabel, id);
				if (existspsd)
					ctfInfo.psd = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				if (existsctf)
					ctfInfo.ctf = md.getValueString(MDLabel.MDL_CTF_MODEL, id);
				if (hasDefocusU)
				    ctfInfo.defocusU = new Float(md.getValueDouble(MDLabel.MDL_CTF_DEFOCUSV, id));

				micrograph = new SupervisedPickerMicrograph(filename, ctfInfo);

				micrographs.add(micrograph);
			}
			if (micrographs.isEmpty())
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", getMicrographsSelFile()));
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

    }
	
    public void loadEmptyMicrographsFromSqlite()
    {
    	if (micrographs == null)
			micrographs = new ArrayList<SupervisedPickerMicrograph>();
		else
			micrographs.clear();
		SupervisedPickerMicrograph micrograph;
		String  filename;
		try
		{
			String selfile = getMicrographsSelFile();
			ScipionMetaData md = new ScipionMetaData(selfile);

			ColumnInfo ci;
			if(md.getSelf().equals("Micrograph"))
			{
				ci = md.getColumnInfo("_filename");
				long[] ids = md.findObjects();
				for (long id : ids)
				{

					filename = md.getValueString(ci.label, id);
					micrograph = new SupervisedPickerMicrograph(filename, null);
					micrographs.add(micrograph);
				}
			}
			else
				throw new IllegalArgumentException(String.format("Labels MDL_MICROGRAPH or MDL_IMAGE not found in metadata %s", selfile));
			
			if (micrographs.isEmpty())
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", getMicrographsSelFile()));
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
    }
	
	
	@Override
	public List<SupervisedPickerMicrograph> getMicrographs()
	{
		return micrographs;
	}

	public synchronized void saveData(Micrograph m)
	{

        JMetaDataIO.saveData((SupervisedPickerMicrograph) m, getOutputPath(m.getPosFile()));
        saveTemplates();

	}

	public SupervisedPickerMicrograph getMicrograph()
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
		this.micrograph = (SupervisedPickerMicrograph) m;
	}

	@Override
	public boolean isValidSize(JFrame parent, int size)
	{


        if (size > ParticlePicker.sizemax) {
             XmippDialog.showInfo(parent, String.format("Max size is %s, %s not allowed", ParticlePicker.sizemax, size));
             return false;
        }
        String outmsg = "";
        int count;
        boolean valid = true;
        List<ManualParticle> outlist = new ArrayList<ManualParticle>();
        for (SupervisedPickerMicrograph m : micrographs)
        {
            count = 0;
            for (ManualParticle p : m.getParticles())
                    if (!m.fits(p.getX(), p.getY(), size))
                    {
                            outlist.add(p);
                            count ++;
                            valid = false;
                    }
            if(count > 0)
                outmsg += String.format("%s %s from micrograph %s will be dismissed if size %s is set.", count, (count == 1)? "particle": "particles", m.getName(), size);
        }
        if(valid)
            return true;
        else
        {
            String msg = outmsg + " Are you sure you want to continue?";
            valid = XmippDialog.showQuestion(parent, msg);
            if(valid)
                for(ManualParticle p: outlist)
                    ((SupervisedPickerMicrograph)p.getMicrograph()).removeParticle(p, this);
            return valid;
                    
        }
	}

        @Override
	public void loadConfig()
	{
		super.loadConfig();
		String file = configfile;
		if (!new File(file).exists())
			return;

		try
		{
			MetaData md = new MetaData(file);
			Mode configmode;
			boolean hasautopercent = md.containsLabel(MDLabel.MDL_PICKING_AUTOPICKPERCENT);
            long id = md.firstObject();

            if(hasautopercent) 
                autopickpercent = md.getValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, id);
            
            dtemplatesnum = md.getValueInt(MDLabel.MDL_PICKING_TEMPLATES, id);
            if (dtemplatesnum == 0)
                    dtemplatesnum = 1;// for compatibility with previous projects
            configmode = Mode.valueOf(md.getValueString(MDLabel.MDL_PICKING_STATE, id));
            isautopick = configmode == Mode.Supervised || configmode == Mode.Review;
            
            if (mode == Mode.Review && configmode == Mode.Manual)//Review mode makes no sense if manual mode
                throw new IllegalArgumentException("Cannot review picking in manual mode, use manual mode instead");

            if (mode != Mode.ReadOnly && mode != Mode.Review)
            	mode = configmode;
            
			md.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}
        
    public boolean isAutopick()
    {
        return isautopick;
    }
    
    @Override
    public synchronized void saveConfig()
    {
        JMetaDataIO.saveConfig(this, configfile);
    }
        
        @Override
	protected synchronized void saveConfig(MetaData md, long id)
	{
		try
		{
			super.saveConfig(md, id);
			md.setValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, getAutopickpercent(), id);
			md.setValueInt(MDLabel.MDL_PICKING_TEMPLATES, getTemplatesNumber(), id);
			md.setValueString(MDLabel.MDL_PICKING_STATE, mode.toString(), id);
                        md.setValueInt(MDLabel.MDL_PICKING_MANUALPARTICLES_SIZE, getManualParticlesNumber(), id);
                        md.setValueInt(MDLabel.MDL_PICKING_AUTOPARTICLES_SIZE, getAutomaticParticlesNumber(), id);

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
                if(autopickpercent == 0)
                    autopickpercent = 50;
		return autopickpercent;
	}

	public void setMode(Mode mode)
	{
		if (mode == Mode.Supervised && getManualParticlesNumber() < classifier.getTrainingParticlesMinimum())
			throw new IllegalArgumentException(String.format("You should have at least %s particles to go to %s mode", classifier.getTrainingParticlesMinimum(), Mode.Supervised));
		this.mode = mode;

		// If supervised
		if(mode == Mode.Supervised)
            isautopick = true;
		if (mode == Mode.Manual)
			convertAutomaticToManual();

	}

	protected void convertAutomaticToManual()
	{
		ManualParticle p;
		for (SupervisedPickerMicrograph m : micrographs)
		{
			m.deleteBelowThreshold();
			for (AutomaticParticle ap : m.getAutomaticParticles())
			{
				if (!ap.isDeleted())
				{
					p = new ManualParticle(ap.getX(), ap.getY(), this, m);
					m.addManualParticle(p, this);
				}
			}

			m.getAutomaticParticles().clear();
			if (m.hasManualParticles())
				m.setState(MicrographState.Manual);
			else
				m.setState(MicrographState.Available);
		}

		saveAllData();
		updateTemplatesStack(false);
                
	}

	public int getManualParticlesNumber()
	{
		int count = 0;
		for (SupervisedPickerMicrograph m : micrographs)
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
			for (SupervisedPickerMicrograph m : micrographs)
			{
				if (m.isEmpty())
					continue;
				for (ManualParticle p : m.getParticles())
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

            this.setSaved(false);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public int getAutomaticParticlesNumber(double threshold)
	{
		int count = 0;
		for (SupervisedPickerMicrograph m : micrographs)
			count += m.getAutomaticParticlesNumber(threshold);
		return count;
	}
        
    public int getAutomaticParticlesNumber()
	{
		int count = 0;
		for (SupervisedPickerMicrograph m : micrographs)
			count += m.getAutomaticParticlesNumber(m.getThreshold());
		return count;
	}

	public void loadAutomaticParticles(SupervisedPickerMicrograph m, MetaData md)
	{

		int x, y;
		AutomaticParticle particle;
		Double cost;
		boolean deleted;
		try
		{
			boolean added = false;

			for (long id : md.findObjects())
			{

				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				cost = md.getValueDouble(MDLabel.MDL_COST, id);
				if (cost == null)
					throw new IllegalArgumentException("Invalid format for " + md.getFilename());
				deleted = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? false : true;
				particle = new AutomaticParticle(x, y, this, m, cost, deleted);
				m.addAutomaticParticle(particle);
				added = true;
			}
			// Change the state of the micrograph if some automatic particles
			// were added
			if (added && m.getState() == MicrographState.Available)
				m.setState(MicrographState.Auto);

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void loadAutomaticParticles(SupervisedPickerMicrograph m, String file)
	{
		if (new File(file).exists() && MetaData.containsBlock(file, particlesAutoBlock))// todo:																						// exists
		{
			MetaData md = new MetaData(getParticlesAutoBlock(file));
			loadAutomaticParticles(m, md);
			md.destroy();
		}
	}

	/** Return the number of particles imported */
	public String importParticlesFromFolder(String path, Format f, String preffix, String suffix, float scale, boolean invertx, boolean inverty)
	{

		if (f == Format.Auto)
		{
			StringBuffer suffixPtr = new StringBuffer(suffix);
			f = detectFormat(path, preffix, suffixPtr);
			suffix = suffixPtr.toString();
		}
		if (f == Format.None)
			throw new IllegalArgumentException("Unable to detect format. You may try specifying format or renaming files");

		String result = "";

		String particlesfile = getExportFile(path);
		if (particlesfile != null && f == Format.Xmipp30)
		{
			importAllParticles(particlesfile);
			return "";
		}

		importSize(path, f, scale);
                String file;
                
		for (SupervisedPickerMicrograph m : micrographs)
		{
            resetMicrograph(m);
            file = Filename.join(path, preffix + m.getName() + suffix);
            if(new File(file).exists())
                result += importParticlesFromFile(file, f, m, scale, invertx, inverty);
		}
                
                
                        
		return result;
	}// function importParticlesFromFolder
        
        

	/** Return the number of particles imported from a file */
	public String importParticlesFromFile(String path, Format f, SupervisedPickerMicrograph m, float scale, boolean invertx, boolean inverty)
	{
		String result = "";
		try
		{
			MetaData md = new MetaData();
			fillParticlesMdFromFile(path, f, m, md, scale, invertx, inverty);
			if (md.size() > 0)
			{
				//Be sure that width and height are loaded
				result = importParticlesFromMd(m, md);
			}
			md.destroy();
		}
		catch (Exception ex)
		{
            XmippDialog.showInfo(null, ex.getMessage());
		}
		return result;
	}// function importParticlesFromFile

	public String importParticlesFromMd(SupervisedPickerMicrograph m, MetaData md)
	{
		long[] ids = md.findObjects();
		int x, y;
		String result = "";
		int size = getSize();

		for (long id : ids)
		{
			x = md.getValueInt(MDLabel.MDL_XCOOR, id);
			y = md.getValueInt(MDLabel.MDL_YCOOR, id);
			if ( !m.fits(x, y, size))// ignore out of bounds particles
			{
				result += XmippMessage.getOutOfBoundsMsg(String.format("Particle on x:%s y:%s\n", x, y)) + " dismissed\n";
				continue;
			}
			m.addManualParticle(new ManualParticle(x, y, this, m, 2), this);

		}
		saveData(m);
                
		return result;
	}// function importParticlesFromMd

	public boolean isReviewFile(String file)
	{
		if (!file.endsWith("xmd"))
			return false;
		boolean result = false;
		boolean nullMicrograph = false;
		try
		{
			MetaData md = new MetaData(file);

			if (md.containsLabel(MDLabel.MDL_XCOOR) && md.containsLabel(MDLabel.MDL_YCOOR) && md.containsLabel(MDLabel.MDL_COST))
			{
				String[] blocksArray = MetaData.getBlocksInMetaDataFile(file);
				List<String> blocks = Arrays.asList(blocksArray);

				String name;
				for (String block : blocks)
				{
					name = block.replace("mic_", "");
					if (getMicrograph(name) == null)
					{
						nullMicrograph = true;
						break;
					}
				}
			}
			md.destroy();

			if (!nullMicrograph)
				result = true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return result;
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

			for (SupervisedPickerMicrograph m : micrographs)
			{
				resetMicrograph(m);
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
			throw new IllegalArgumentException(e.getMessage());
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

			for (SupervisedPickerMicrograph m : micrographs)
			{
				resetMicrograph(m);
				block = "mic_" + m.getName();
				if (blocks.contains(block))
				{
					String blockName = block + "@" + file;
					md.read(blockName);
					width = (int) (m.getWidth() / scale);// original width
					height = (int) (m.getHeigth() / scale);// original height
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
			throw new IllegalArgumentException(e.getMessage());
		}
		return null;

	}// function importAllParticles



	public boolean hasManualParticles()
	{
		for (SupervisedPickerMicrograph m : micrographs)
		{
			if (m.hasManualParticles())
				return true;
		}
		return false;
	}

	public synchronized int getTemplateIndex()
	{
		return templateindex;
	}


	public void loadMicrographData(SupervisedPickerMicrograph m)
	{
		try
		{
			String file = getOutputPath(m.getPosFile());
			MicrographState state;
			int micautopickpercent;
			if (!new File(file).exists())
				return;
			if (MetaData.containsBlock(file, "header"))
			{
				MetaData md = new MetaData("header@" + file);
				boolean hasautopercent = md.containsLabel(MDLabel.MDL_PICKING_AUTOPICKPERCENT);
				long id = md.firstObject();
				state = MicrographState.valueOf(md.getValueString(MDLabel.MDL_PICKING_MICROGRAPH_STATE, id));
				micautopickpercent = hasautopercent ? md.getValueInt(MDLabel.MDL_PICKING_AUTOPICKPERCENT, id) : getAutopickpercent();
				double threshold = md.getValueDouble(MDLabel.MDL_COST, id);
				m.setThreshold(threshold);
				m.setState(state);
				m.setAutopickpercent(micautopickpercent);
				md.destroy();
			}
			loadManualParticles(m, file);
			loadAutomaticParticles(m, file);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void loadManualParticles(SupervisedPickerMicrograph m, String file)
	{
		if (new File(file).exists() && MetaData.containsBlock(file, getParticlesBlockName(Format.Xmipp301)))
		{

			int x, y;
			ManualParticle particle;

			try
			{
				MetaData md = new MetaData(getParticlesBlock(file));

				for (long id : md.findObjects())
				{
					x = md.getValueInt(MDLabel.MDL_XCOOR, id);
					y = md.getValueInt(MDLabel.MDL_YCOOR, id);
					particle = new ManualParticle(x, y, this, m);
					m.addManualParticle(particle, this);
				}
				md.destroy();
			}
			catch (Exception e)
			{
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
	 * @param trainmic
	 */
	public void trainAndAutopick(SupervisedPickerJFrame frame, SupervisedPickerMicrograph trainmic)
	{
		frame.getCanvas().setEnabled(false);
		XmippWindowUtil.blockGUI(frame, "Training and Autopicking...");
		// Change the state of the micrograph and save changes
		micrograph.setState(MicrographState.Supervised);
		saveData(micrograph);
        setChanged(false);
		// Create the metadata with each micrograph and its .pos file
		// Add all micrographs with manual particles except the train
		// micrograph
        ArrayList<MDRow> trainRows = new ArrayList<MDRow>();
		for (SupervisedPickerMicrograph m : micrographs)
		{
			if (m.hasManualParticles() && !m.equals(trainmic))
				addMicrographPos(trainRows, m);
		}
		// Add the train micrograph as the last one
			
		addMicrographPos(trainRows, trainmic);
                
		new Thread(new TrainRunnable(frame, trainRows.toArray(new MDRow[]{}), trainmic), "TrainAndAutopick").start();
	}

	/** Helper function to add a micrograph and its .pos file */
	private void addMicrographPos(List<MDRow> trainRows, SupervisedPickerMicrograph micrograph)
	{
		String posfile = getOutputPath(micrograph.getPosFile());
                MDRow row = new MDRow();
		row.setValueString(MDLabel.MDL_MICROGRAPH, micrograph.getFile());
		row.setValueString(MDLabel.MDL_MICROGRAPH_PARTICLES, posfile);
                trainRows.add(row);

	}

    @Override
    public int getParticlesCount() {
        return getManualParticlesNumber() + getAutomaticParticlesNumber();
    }

    public void loadParticles(SupervisedPickerMicrograph micrograph, Particle[] autopickRows)
    {
            Rectangle rectangle = micrograph.getRectangle();
            int x, y;
            double cost;
            AutomaticParticle ap;
            for (Particle row: autopickRows)
            {
                    x = row.getX();
                    y = row.getY();
                    cost = row.getCost();
                    ap = new AutomaticParticle(x, y, SupervisedParticlePicker.this, micrograph, cost, false);
                    if (rectangle != null &&  rectangle.contains(new Point(x, y)))
                    {
                            ap.setDeleted(true);
                            ap.setCost(-1);
                    }
                    micrograph.addAutomaticParticle(ap);
            }
    }


	public void autopick(SupervisedPickerJFrame frame, SupervisedPickerMicrograph next)
	{
                frame.getCanvas().setEnabled(false);
		XmippWindowUtil.blockGUI(frame, "Autopicking...");

        try {
            new Thread(new AutopickRunnable(frame, next), "Autopick").start();
        }catch (Exception e){
            this.logger.log(Level.SEVERE,"Autopick thread failed.");
        }


	}

    public Classifier getClassifier() {
        return classifier;
    }

	

	public void correctAndAutopick(SupervisedPickerJFrame frame, SupervisedPickerMicrograph current, SupervisedPickerMicrograph next)
	{
		current.setState(MicrographState.Corrected);
		if (getMode() == Mode.Supervised && next.getState() == MicrographState.Available)
                    next.setState(MicrographState.Supervised);
		saveData(current);


        MDRow[] autoRows = getAutomaticRows(current);
        MDRow[] addedRows = getAddedRows(current);
		frame.getCanvas().setEnabled(false);

        boolean isautopick = getMode() == Mode.Supervised && next.getState() == MicrographState.Available;
        // Show progress wheel
        if(isautopick)
            XmippWindowUtil.blockGUI(frame, "Correcting and Autopicking...");
        else
            XmippWindowUtil.blockGUI(frame, "Correcting ...");

		try {
            new Thread(new CorrectAndAutopickRunnable(frame, addedRows, autoRows, next), "CorrectAndAutopick").start();
        }catch (Exception e){
		    this.logger.log(Level.SEVERE,"Correct and Autopick thread failed.");
        }

	}

	private MDRow[] getAddedRows(SupervisedPickerMicrograph m)
	{
                Rectangle correctout = m.getRectangle();
                ArrayList<MDRow> rows = new ArrayList<MDRow>();
                MDRow row;
                for (ManualParticle p : m.getManualParticles())
                        if (correctout == null || !correctout.contains(new Point(p.getX(), p.getY())))
                        {
                                row = new MDRow();
                                row.setValueInt(MDLabel.MDL_XCOOR, p.getX());
                                row.setValueInt(MDLabel.MDL_YCOOR, p.getY());
                                rows.add(row);
                        }
		return rows.toArray(new MDRow[]{});
	}
    // Returns all automatic rows that matches the conditions, threshold is taken into account.
    private MDRow[] getAutomaticRows(SupervisedPickerMicrograph m)
	{
                ArrayList<MDRow> rows = new ArrayList<MDRow>();
                MDRow row;
		try
		{

			for (AutomaticParticle p : m.getAutomaticParticles())
			{
			    row = new MDRow();
				row.setValueInt(MDLabel.MDL_XCOOR, p.getX());
				row.setValueInt(MDLabel.MDL_YCOOR, p.getY());
				row.setValueDouble(MDLabel.MDL_COST, p.getCost());
                // row.setValueInt(MDLabel.MDL_ENABLED, p.isDeleted()? -1: 1);
                row.setValueInt(MDLabel.MDL_ENABLED, p.isUnavailable()? -1: 1);
				rows.add(row);
			}
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
		return rows.toArray(new MDRow[]{});
	}

	
	public void resetMicrograph(SupervisedPickerMicrograph micrograph)
	{
		micrograph.getManualParticles().clear();
		micrograph.getAutomaticParticles().clear();
		micrograph.setState(MicrographState.Available);
		micrograph.setThreshold(0);
		micrograph.resetParticlesRectangle();
		new File(getOutputPath(micrograph.getPosFile())).delete();

	}

	//This method will only be called from the interface if xmipp picker is used
	public void correct()
	{
		micrograph.setState(MicrographState.Corrected);
		saveData(micrograph);
		MDRow[] manualRows = getAddedRows(micrograph);
		MDRow[] autoRows = getAutomaticRows(micrograph);
		((PickingClassifier)classifier).correct(manualRows, autoRows);
	}
        
    public int getParticlesThreshold() {
        return classifier.getTrainingParticlesMinimum();
    }
    
     
    public void updateTemplatesStack(boolean resize)
    {
        updateTemplatesStack(getTemplatesNumber(), resize);
    }
     
    public void updateTemplatesStack(int num, boolean resize)
    {
        if(mode != Mode.Manual)
            return;
        initTemplates(num);//to avoid any error with templates while updating
        new UpdateTemplatesTask(resize).execute();
    }
    
    public void setTemplatesDialog(TemplatesJDialog d)
    {
            dialog = d;
    }
    
    public class UpdateTemplatesTask extends SwingWorker<String, Object>
    {
        private final boolean resize;
        public UpdateTemplatesTask(boolean resize)
        {
            this.resize = resize;
        }

        @Override
        protected String doInBackground() throws Exception
        {
            try
            {
                if(uttask != null)
                    uttask.cancel(true);
                uttask = this;


                ImageGeneric templates = new ImageGeneric(ImageGeneric.Float);
                templates.resize(getSize(), getSize(), 1, getTemplatesNumber());

                ImageGeneric igp;
                List<ManualParticle> particles;

                ManualParticle particle;
                double[] align;

                //FIXME: the template update needs to be done in a
                // more efficient way, now we are using this maxcount
                // to limit the number of particles used in the update
                int count = 0;
                int maxcount = 50;
                //FIXME: This is a tweak to avoid losing time with big particles
                if (getSize() > 256)
                        maxcount = 10;
                float[] matrix;
                int templateindex = 0;
                for (SupervisedPickerMicrograph m : getMicrographs())
                {
                        particles = m.getManualParticles();
                        for (int i = 0; i < particles.size(); i++)
                        {
                                if(isCancelled())
                                    return null;
                                if(count < maxcount)
                                {
                                    particle = particles.get(i);
                                    if(resize)
                                        particle.resetImage();
                                    igp = particle.getImageGeneric();
                                    if (templateindex < getTemplatesNumber())
                                    {
                                            matrix = igp.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.FIRST_SLICE);
                                            templates.setArrayFloat(matrix, ImageGeneric.FIRST_IMAGE + templateindex, ImageGeneric.FIRST_SLICE);
                                            templateindex ++;
                                    }
                                    else
                                    {
                                            align = templates.alignImage(igp);
                                            particle.setLastalign(align);
                                            templates.applyAlignment(igp, particle.getTemplateIndex(), particle.getTemplateRotation(), particle.getTemplateTilt(), particle.getTemplatePsi());
                                    }
                                    count++;
                                }

                        }
                }

                updateTemplates(templates, templateindex);
            }
            catch (Exception e)
            {
                    e.printStackTrace();
                    throw new IllegalArgumentException(e);
            }
            return "";
        }

        @Override
        protected void done() {
                    if(dialog != null && dialog.isVisible())
                                   dialog.loadTemplates(true);
        }

    }
    
    public boolean isCorrectPending()//TO BE COMPLETED
    {
        return getMode() == Mode.Supervised && getMicrograph().getState() == MicrographState.Supervised && isChanged();
    }

	public Color getAutomaticColor() {
		return getMode() == Mode.Automatic ? getColor() : moveColorHue(getColor(), 0.66f);
	}

	public Color getDeletedColor() {

		return moveColorHue(getColor(), 0.33f);
	}

	private Color moveColorHue(Color color, float hueValue){
		// Get saturation and brightness.
		float[] hsbVals = new float[3];
		Color.RGBtoHSB(color.getRed(), color.getGreen(), color.getBlue(), hsbVals);

		// Pass .5 (= 180 degrees) as HUE
		float newHue = hsbVals[0] - hueValue;

		if (newHue < 0) newHue = newHue + 1f;

		Color newColor = new Color(Color.HSBtoRGB(newHue, hsbVals[1], hsbVals[2]));
		return newColor;
	}

}
