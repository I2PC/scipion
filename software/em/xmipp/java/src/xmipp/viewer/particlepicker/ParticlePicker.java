package xmipp.viewer.particlepicker;

import ij.plugin.frame.Recorder;
import java.awt.Color;
import java.io.File;
import java.io.FilenameFilter;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.jni.Program;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.training.model.Mode;

/**
 * Business object for ParticlePicker common GUI. SingleParticlePicker and
 * TiltPairPicker inherit from this class.
 *
 * @author airen
 *
 */
public abstract class ParticlePicker {

    protected String macrosfile;
    protected static Logger logger;
    protected String outputdir = ".";
    protected boolean changed;
    protected Mode mode;
    protected List<IJCommand> filters;
    protected String selfile;
    
    protected String configfile;
    protected final String[] availableFilters = new String[]{"Duplicate", "Bandpass Filter...", "Anisotropic Diffusion...", "Mean Shift",
                   "Subtract Background...", "Gaussian Blur...", "Brightness/Contrast...", "Invert LUT"};

    
    static String xmippsmoothfilter = "Xmipp Smooth Filter";
    public static final String particlesAutoBlock = "particles_auto";

    protected Color color;
    protected int size;

    public static final int sizemax = 2000;
    protected String block;

    protected static Color[] colors = new Color[]{Color.BLUE, Color.CYAN, Color.GREEN, Color.MAGENTA, Color.ORANGE, Color.PINK, Color.YELLOW};

    protected static int nextcolor;
    public final HashMap<Format, String> emextensions;
    protected ParticlePickerParams params;
    protected static ParticlePicker picker;

    public static Color getNextColor() {
        Color next = colors[nextcolor];
        nextcolor++;
        if (nextcolor == colors.length) {
            nextcolor = 0;
        }
        return next;
    }

    public static String getParticlesBlock(Format f, String file) {
        String blockname = getParticlesBlockName(f);
        if (f == null) {
            return null;
        }
        return blockname + "@" + file;

    }

    public static String getParticlesBlockName(Format f) {
        if (f == Format.Xmipp301) {
            return "particles";
        }
        if (f == Format.Xmipp30) {
            return "DefaultFamily";
        }
        return null;
    }

    public static String getParticlesBlock(String file) {
        return getParticlesBlock(Format.Xmipp301, file);

    }
    
    public ParticlePickerParams getParams()
    {
        return params;
    }
    

    public String getParticlesAutoBlock(String file) {
        return particlesAutoBlock + "@" + file;
    }

    public String getParticlesBlock(Micrograph m) {
        return getParticlesBlock(getOutputPath(m.getPosFile()));

    }

    public String getParticlesAutoBlock(Micrograph m) {
        return getParticlesAutoBlock(getOutputPath(m.getPosFile()));
    }

    public ParticlePicker(String selfile, Mode mode, ParticlePickerParams params) {
        this(null, selfile, ".", mode, params);
    }

    public ParticlePicker(String selfile, String outputdir, Mode mode, ParticlePickerParams params) {
        this(null, selfile, outputdir, mode, params);
    }

    public ParticlePicker(String block, String selfile, String outputdir, Mode mode, ParticlePickerParams params) {
    	
        this.params = params;
        this.block = block;
        this.outputdir = outputdir;
        if(!new File(outputdir).isDirectory())
            throw new IllegalArgumentException(XmippMessage.getInvalidDirectoryMsg(outputdir));
        this.selfile = selfile;
        this.mode = mode;
        this.configfile = getOutputPath("config.xmd");
        picker = this;
        filters = new ArrayList<IJCommand>();
        loadEmptyMicrographs();
        loadConfig();
        emextensions = new HashMap<Format, String>();
        emextensions.put(Format.Xmipp, ".pos");
        emextensions.put(Format.Xmipp24, ".pos");
        emextensions.put(Format.Xmipp30, ".pos");
        emextensions.put(Format.Xmipp301, ".pos");
        emextensions.put(Format.Relion, ".star");
        emextensions.put(Format.Eman, ".box");
        
    }

    public void loadConfig() {

        String file = configfile;
        filters.clear();
        if (!new File(file).exists()) {
            color = getNextColor();
            size = getDefaultSize();
            setMicrograph(getMicrographs().get(0));
            filters.add(new IJCommand("Gaussian Blur...", "sigma=2"));
            saveConfig();
            return;
        }

        String mname;
        int rgb;
        try {
            MetaData md = new MetaData("properties@" + file);
            for (long id : md.findObjects()) {

                mname = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
                setMicrograph(getMicrograph(mname));
                rgb = md.getValueInt(MDLabel.MDL_COLOR, id);
                color = new Color(rgb);
                size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
            }
            md.destroy();
            String command, options;
            if (MetaData.containsBlock(file, "filters")) {
                md = new MetaData("filters@" + file);
                long[] ids = md.findObjects();
                for (long id : ids) {
                    command = IJCommand.getString(md.getValueString(MDLabel.MDL_MACRO_CMD, id));
                    options = IJCommand.getString(md.getValueString(MDLabel.MDL_MACRO_CMD_ARGS, id));
                    if (options.equals("NULL")) {
                        options = "";
                    }
                    filters.add(new IJCommand(command, options));
                }
                md.destroy();
            }
        } catch (Exception e) {
            getLogger().log(Level.SEVERE, e.getMessage(), e);
            throw new IllegalArgumentException(e.getMessage());
        }

    }

    public synchronized void saveConfig() {
        try {
            MetaData md;
            String file = configfile;
            md = new MetaData();
            long id = md.addObject();
            saveConfig(md, id);
            md.write("properties@" + file);
            md.destroy();
            String cmd, options;
            md = new MetaData();
            for (IJCommand f : filters) {
                id = md.addObject();
                cmd = f.getMdCommand();
                md.setValueString(MDLabel.MDL_MACRO_CMD, cmd, id);
                options = f.getMdOptions();
                md.setValueString(MDLabel.MDL_MACRO_CMD_ARGS, options, id);
            }
            md.writeBlock("filters@" + file);
            md.destroy();
            

        } catch (Exception e) {
            getLogger().log(Level.SEVERE, e.getMessage(), e);
            throw new IllegalArgumentException(e.getMessage());
        }

    }

    protected void saveConfig(MetaData md, long id) {
        md.setValueString(MDLabel.MDL_MICROGRAPH, getMicrograph().getName(), id);
        md.setValueInt(MDLabel.MDL_COLOR, getColor().getRGB(), id);
        md.setValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, getSize(), id);
    }
    
    public Format detectFileFormat(String path) {
        if(!new File(path).exists())
            return Format.None;
        if (path.endsWith(emextensions.get(Format.Xmipp)))
        {
            if(MetaData.isPlainPos(path) ) 
                return Format.Xmipp24;
            if (MetaData.containsBlock(path, getParticlesBlockName(Format.Xmipp30))) {
                return Format.Xmipp30;
            }
            return Format.Xmipp301;
        }
        else if (path.endsWith(emextensions.get(Format.Eman))) {
            return Format.Eman;
        }
        else if (path.endsWith(emextensions.get(Format.Relion)))
            return Format.Relion;
         
        return Format.None;
    }

    public Format detectFormat(String path, String preffix, String suffix) {
        String file;
        Format f;
        for (Micrograph m : getMicrographs())
        {
            file = Filename.join(path, preffix + m.getName() + suffix);
            f = detectFileFormat(file);
            if(f != Format.None)
                return f;
        }
        return Format.None;
    }

    public String getExportFile(String path) {
        String particlesfile;
        File folderToScan = new File(path);
        File[] listOfFiles = folderToScan.listFiles();

        //Checking file with exported particles
        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile()) {
                particlesfile = listOfFiles[i].getName();
                if (particlesfile.endsWith("extract_list.xmd")) {
                    return Filename.join(path, particlesfile);
                }
            }
        }
        return null;
    }


    public String getPosFileFromXmipp24Project(String projectdir, String mname) {
        String suffix = ".raw.Common.pos";
        return String.format("%1$s%2$sPreprocessing%2$s%3$s%2$s%3$s%4$s", projectdir, File.separator, mname, suffix);
    }

    public abstract void loadEmptyMicrographs();



    protected boolean isRegisteredFilter(String command2) {
        for (IJCommand f : filters) {
            if (f.getCommand().equals(command2)) {
                return true;
            }
        }
        return false;
    }

    public void updateFilters(String command) {
        if (command != null) {
            String options = "";
            if (Recorder.getCommandOptions() != null) {
                options = Recorder.getCommandOptions();
            }

            if (!isFilterSelected(command) && Arrays.asList(availableFilters).contains(command)) {
                addFilter(command, options);
            } else if (!(options == null || options.equals(""))) {
                for (IJCommand f : filters) {
                    if (f.getCommand().equals(command)) {
                        f.setOptions(options);
                    }
                }
            }

            saveConfig();
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

    public Mode getMode() {
        return mode;
    }

    public static Logger getLogger() {
        try {
            if (logger == null) {
				//				FileHandler fh = new FileHandler("PPicker.log", true);
                //				fh.setFormatter(new SimpleFormatter());
                logger = Logger.getLogger("PPickerLogger");
                // logger.addHandler(fh);
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

    public abstract List<? extends Micrograph> getMicrographs();

    public int getMicrographIndex() {
        return getMicrographs().indexOf(getMicrograph());
    }

    public String getTemplatesFile(String name) {
        return getOutputPath(name + "_templates.stk");
    }

    public Mode validateState(Mode state) {

        if (mode == Mode.Review && state != Mode.Review) {
            setChanged(true);
            return Mode.Review;
        }
        if (mode == Mode.Manual && !(state == Mode.Manual || state == Mode.Available)) {
            throw new IllegalArgumentException(String.format("Can not use %s mode on this data", mode));
        }
        if (mode == Mode.Supervised && state == Mode.Review) {
            throw new IllegalArgumentException(String.format("Can not use %s mode on this data", mode));
        }
        return state;

    }// function validateState

    public void saveData() {
        saveConfig();
		//saveData(getMicrograph());
        //setChanged(false);
    }// function saveData

    public abstract void saveData(Micrograph m);

    public Micrograph getMicrograph(String name) {
        for (Micrograph m : getMicrographs()) {
            if (m.getName().equalsIgnoreCase(name)) {
                return m;
            }
        }
        return null;
    }

    void removeFilter(String filter) {
        for (IJCommand f : filters) {
            if (f.getCommand().equals(filter)) {
                filters.remove(f);
                saveConfig();
                break;
            }
        }
    }// function removeFilter

    public boolean isFilterSelected(String filter) {
        for (IJCommand f : filters) {
            if (f.getCommand().equals(filter)) {
                return true;
            }
        }
        return false;
    }// function isFilterSelected

    

    /**
     * Return the number of particles imported from a file
     */
    public void fillParticlesMdFromFile(String path, Format f, Micrograph m, MetaData md, float scale, boolean invertx, boolean inverty) {
        
        if (f == Format.Auto) {
            f = detectFileFormat(path);
        }
        String blockname;
        switch (f) {
            case Xmipp24:
                md.readPlain(path, "xcoor ycoor");
                break;
            case Xmipp30:
                blockname = getParticlesBlockName(f); //only used with Xmipp
                if (MetaData.containsBlock(path, blockname)) {
                    md.read(getParticlesBlock(f, path));
                }
                break;
            case Xmipp301:
                fillParticlesMdFromXmipp301File(path, m, md);
                break;
            case Eman:
                fillParticlesMdFromEmanFile(path, m, md, scale);
                break;
             case Relion:
                fillParticlesMdFromRelionFile(path, m, md);
                
                break;
            default:
                md.clear();
        }
        if (scale != 1.f) {
            String command = String.format(Locale.ENGLISH, "xcoor=xcoor*%f,ycoor=ycoor*%f", scale, scale);
            md.operate(command);
        } 
        int width = (int) (m.getWidth() / scale);// original width
        int height = (int) (m.getHeigth() / scale);// original height
        if (invertx) {
            md.operate(String.format("xcoor=%d-xcoor", width));
        }
        if (inverty) {
            md.operate(String.format("ycoor=%d-ycoor", height));
        }
       
        
    }// function importParticlesFromFile
    
    public void fillParticlesMdFromXmipp301File(String file, Micrograph m, MetaData md) {
        
         String blockname = getParticlesBlockName(Format.Xmipp301); //only used with Xmipp
                if (MetaData.containsBlock(file, blockname)) {
                    md.read(getParticlesBlock(Format.Xmipp301, file));
                }
                if (MetaData.containsBlock(file, particlesAutoBlock)) {
                    MetaData mdAuto = new MetaData();
                    mdAuto.read(getParticlesAutoBlock(file));
                    mdAuto.removeDisabled();
                    md.unionAll(mdAuto);
                    mdAuto.destroy();
                }
                
    }

    public void fillParticlesMdFromRelionFile(String file, Micrograph m, MetaData md) {
         MetaData relionmd = new MetaData(file);
         long[] ids = relionmd.findObjects();
         int xcoor, ycoor;
         long id;
         for(int i = 0; i < ids.length; i++ )
         {
            xcoor = (int)relionmd.getValueDouble(MDLabel.RLN_IMAGE_COORD_X, ids[i]);
            ycoor = (int)relionmd.getValueDouble(MDLabel.RLN_IMAGE_COORD_Y, ids[i]);
            id = md.addObject();
            md.setValueInt(MDLabel.MDL_XCOOR, xcoor, id);
            md.setValueInt(MDLabel.MDL_YCOOR, ycoor, id);
         }
         relionmd.destroy();
                
    }
    
    public void fillParticlesMdFromEmanFile(String file, Micrograph m, MetaData md, float scale) {
        // inverty = true;
        md.readPlain(file, "xcoor ycoor particleSize");
		//System.err.format("After readPlain: md.size: %d\n", md.size());

        long fid = md.firstObject();
        int size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, fid);

        int half = size / 2;
        //System.err.format("Operate string: %s\n", String.format("xcoor=xcoor+%d,ycoor=ycoor+%d", half, half));
        md.operate(String.format("xcoor=xcoor+%d,ycoor=ycoor+%d", half, half));

        setSize(Math.round(size * scale));
    }// function fillParticlesMdFromEmanFile

    public void runXmippProgram(String program, String args) {
        try {
            Program.runByName(program, args);
        } catch (Exception e) {
            getLogger().log(Level.SEVERE, e.getMessage(), e);
            throw new IllegalArgumentException(e);
        }
    }

    public abstract Micrograph getMicrograph();

    public abstract void setMicrograph(Micrograph m);

    public abstract boolean isValidSize(JFrame parent, int size);

    public int getSize() {
        return size;
    }

    public void setSize(int size) {
       
        this.size = size;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
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
    
    public void importSize(String path, Format f, double scale)
    {
        if(f == Format.Xmipp301)
                {
                    String configfile = String.format("%s%s%s", path, File.separator, "config.xmd");
                    if(new File(configfile).exists())
                    {
                        MetaData configmd = new MetaData(configfile);
                        int size = configmd.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, configmd.firstObject());
                        size = (int)(size * scale);
                        setSize(size);
                        saveConfig();
                        configmd.destroy();
                    }
                }
    }

    public static File[] getCoordsFiles(String folder)
        {
            return new File(folder).listFiles(new FilenameFilter() {
                    public boolean accept(File dir, String name) {
                        String lowercasename = name.toLowerCase();
                        return lowercasename.endsWith(".pos") || lowercasename.endsWith(".box");
                    }
                });
        }

    
     
     
     
     public boolean isScipionSave()
     {
         if(params == null)
             return false;
         return params.port != null && mode != Mode.ReadOnly;
     }

    public String getScipionSaveCommand() {
        String cmd = String.format("run function registerCoords %s %s", params.dbpath, params.protid);
        return cmd;
    }
     
    
    public abstract int getParticlesCount();

    public String getProtId() {
        return params.protid;
    }

    
    public static ParticlePicker getPicker()
    {
        return picker;
    }
     

   

}
