package xmipp.utils;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import xmipp.jni.ImageGeneric;


/**
 *
 * @author Juanjo Vega
 */
public class Params {

	//Some constants definitions
    public final static String FILE = "i";
    public final static String DIRECTORY = "dir";
    public final static String INPUT_VECTORSFILE = "vectors";
    public final static String INPUT_CLASSESFILE = "classes";
    public final static String INPUT_DATAFILE = "data";
    public final static String POLL = "poll";
    public final static String TABLE_ROWS = "rows";
    public final static String TABLE_COLUMNS = "columns";
    public final static String ZOOM = "zoom";
    public final static String OPENING_MODE = "mode";
    public final static String OPENING_MODE_DEFAULT = "default";
    public final static String OPENING_MODE_IMAGE = "image";
    public final static String OPENING_MODE_GALLERY = "gallery";
    public final static String OPENING_MODE_METADATA = "metadata";
    public final static String RENDER_IMAGES = "render";
    public final static String VISIBLE_LABELS = "visible";
    public final static String ORDER_LABELS = "order";
    public final static String SORTBY_LABEL = "sortby";
    public final static String DISPLAY_LABELS = "labels";
    public final static String FILTER = "filter";
    public final static String FILTERS_SEPARATOR = ",";
    public final static String PORT = "port";
    public final static String DOWNSAMPLING = "downsampling";
    public final static String BAD_PIXELS_FACTOR = "factor";
    public final static String W1 = "w1";
    public final static String W2 = "w2";
    public final static String RAISED_W = "raised_w";
    public final static String MASK = "mask";
    public final static String DEBUG = "debug";
    public final static String MASKTOOLBAR = "mask_toolbar";
    public final static String VIEW = "view";
    public final static String NO_GEO = "dont_apply_geo";
    public final static String NO_WRAP = "dont_wrap";
    public final static String OBJECT_CMDS = "object_commands";
    public final static String SAMPLINGRATE = "sampling_rate";
    public final static String CHIMERAPORT = "chimera_port";
    public final static String INVERTY = "inverty";
    private static final String NO_RECALCULATE_CTF = "dont_recalc_ctf";


    public String directory;
    public String files[];
    public Integer port;
    public boolean singleSelection;
    public String mode = OPENING_MODE_DEFAULT;
    public boolean poll;
    public Integer zoom;
    public String[] renderLabels = new String[]{"first"}; //Label to render, by default first
    public String renderLabel;
    public String[] displayLabels;
    public boolean debug = false;
    public boolean mask_toolbar = false;
    public int rows = -1, columns = -1;
    public double bad_pixels_factor;
    public double w1;
    public double w2;
    public double raised_w;
    public double downsampling;
    public String maskFilename;
    public String vectorsFilename;
    public String classesFilename;
    public String dataFilename;
    public int resliceView = ImageGeneric.Z_NEG; 
    public boolean useGeo = true;
    public boolean wrap = true;
    protected Options options;
    protected CommandLine cmdLine;
    public String[] visibleLabels;
    public String[] orderLabels;
    public String[] sortby;
    private String block;
    public String[] objectCommands;
    private Float samplingRate;
    public Integer chimeraPort;
    public boolean inverty;
	public boolean renderImages = true;
    public boolean recalculateCTF = true;


    public Params() {
    }

    public Params(String args[]) {
        options = new Options();
        defineArgs();
        processArgs(args);
    }

    public void defineArgs()
    {
         options.addOption(FILE, true, "");
        // It should be able to handle multiple files.
        options.getOption(FILE).setOptionalArg(true);
        options.getOption(FILE).setArgs(Integer.MAX_VALUE);

        options.addOption(DIRECTORY, true, "");
        options.addOption(PORT, true, "");

        options.addOption(OPENING_MODE, true, "");
        options.addOption(POLL, false, "");
        options.addOption(ZOOM, true, "");
        Option opt = new Option(RENDER_IMAGES, true, "");
        opt.setArgs(Integer.MAX_VALUE);
        options.addOption(opt);
        opt = new Option(VISIBLE_LABELS, "");
        opt.setArgs(Integer.MAX_VALUE);
        options.addOption(opt);
        opt = new Option(ORDER_LABELS, "");
        opt.setArgs(Integer.MAX_VALUE);
        options.addOption(opt);
        
        opt = new Option(SORTBY_LABEL, "");
        opt.setArgs(2);
        options.addOption(opt);
        
        opt = new Option(DISPLAY_LABELS, "");
        opt.setArgs(Integer.MAX_VALUE);
        options.addOption(opt);
        options.addOption(DEBUG, false, "");
        options.addOption(MASKTOOLBAR, false, "");
        options.addOption(TABLE_ROWS, true, "");
        options.addOption(TABLE_COLUMNS, true, "");

        options.addOption(BAD_PIXELS_FACTOR, true, "");

        options.addOption(W1, true, "");
        options.addOption(W2, true, "");
        options.addOption(RAISED_W, true, "");

        options.addOption(DOWNSAMPLING, true, "");

        options.addOption(MASK, false, "");

        options.addOption(INPUT_VECTORSFILE, true, "");
        options.addOption(INPUT_CLASSESFILE, true, "");
        options.addOption(INPUT_DATAFILE, true, "");
        options.addOption(VIEW, true, "z");
        options.addOption(NO_GEO, false, "");
        options.addOption(NO_WRAP, false, "");
        opt = new Option(OBJECT_CMDS, "");
        opt.setArgs(Integer.MAX_VALUE);
        options.addOption(opt);
        options.addOption(SAMPLINGRATE, true, "");
        options.addOption(CHIMERAPORT, true, "");
        options.addOption(INVERTY, false, "");
        // Do not offer to recalculate CFT if true
        options.addOption(NO_RECALCULATE_CTF, false, "");


    }
    
    
    public void processArgs(String args[]) {


        try {
            BasicParser parser = new BasicParser();
            cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(FILE)) {
                files = cmdLine.getOptionValues(FILE);
            }

            if (cmdLine.hasOption(DIRECTORY)) {
                directory = cmdLine.getOptionValue(DIRECTORY);
            }

            if (cmdLine.hasOption(PORT)) {
                port = Integer.parseInt(cmdLine.getOptionValue(PORT));
            }

            if (cmdLine.hasOption(OPENING_MODE)) {
                mode = cmdLine.getOptionValue(OPENING_MODE);
            } else {
                mode = OPENING_MODE_DEFAULT;
            }

            poll = cmdLine.hasOption(POLL);

            if (cmdLine.hasOption(ZOOM)) {
                zoom = Integer.parseInt(cmdLine.getOptionValue(ZOOM));
            }

            if(cmdLine.hasOption(RENDER_IMAGES))
            {
            	renderLabels = cmdLine.getOptionValues(RENDER_IMAGES);
            	renderImages = !(renderLabels.length == 1 && renderLabels[0].equalsIgnoreCase("no"));
                
            }
            if(cmdLine.hasOption(VISIBLE_LABELS))
            {
            	visibleLabels = cmdLine.getOptionValues(VISIBLE_LABELS);
            }
            if(cmdLine.hasOption(ORDER_LABELS))
            {
            	orderLabels = cmdLine.getOptionValues(ORDER_LABELS);
            }
            if(cmdLine.hasOption(SORTBY_LABEL))
            {
                
            	sortby = cmdLine.getOptionValues(SORTBY_LABEL);
                
            }
            if(cmdLine.hasOption(DISPLAY_LABELS))
            	displayLabels = cmdLine.getOptionValues(DISPLAY_LABELS);
            
            debug = cmdLine.hasOption(DEBUG);
            mask_toolbar = cmdLine.hasOption(MASKTOOLBAR);
            
            if (cmdLine.hasOption(TABLE_ROWS)) {
                rows = Integer.parseInt(cmdLine.getOptionValue(TABLE_ROWS));
            }

            if (cmdLine.hasOption(TABLE_COLUMNS)) {
                columns = Integer.parseInt(cmdLine.getOptionValue(TABLE_COLUMNS));
            }

            if (cmdLine.hasOption(BAD_PIXELS_FACTOR)) {
                bad_pixels_factor = Double.parseDouble(cmdLine.getOptionValue(BAD_PIXELS_FACTOR));
                System.out.println("bad_pixels.-.. "+bad_pixels_factor);
            }

            if (cmdLine.hasOption(W1)) {
                w1 = Double.parseDouble(cmdLine.getOptionValue(W1));
            }

            if (cmdLine.hasOption(W2)) {
                w2 = Double.parseDouble(cmdLine.getOptionValue(W2));
            }

            if (cmdLine.hasOption(RAISED_W)) {
                raised_w = Double.parseDouble(cmdLine.getOptionValue(RAISED_W));
            }

            if (cmdLine.hasOption(DOWNSAMPLING)) {
                downsampling = Double.parseDouble(cmdLine.getOptionValue(DOWNSAMPLING));
            }

            

            if (cmdLine.hasOption(INPUT_VECTORSFILE)) {
                vectorsFilename = cmdLine.getOptionValue(INPUT_VECTORSFILE);
            }

            if (cmdLine.hasOption(INPUT_CLASSESFILE)) {
                classesFilename = cmdLine.getOptionValue(INPUT_CLASSESFILE);
            }

            if (cmdLine.hasOption(INPUT_DATAFILE)) {
                dataFilename = cmdLine.getOptionValue(INPUT_DATAFILE);
            }
            
            
            useGeo = !cmdLine.hasOption(NO_GEO);
            wrap = !cmdLine.hasOption(NO_WRAP);
            recalculateCTF = !cmdLine.hasOption(NO_RECALCULATE_CTF);
            
            if (cmdLine.hasOption(VIEW)){
            	String view = cmdLine.getOptionValue(VIEW);
            	if (view.equals("z"))
            		resliceView = ImageGeneric.Z_NEG;
            	else if (view.equals("y"))
            		resliceView = ImageGeneric.Y_NEG;
            	else if (view.equals("x"))
            		resliceView = ImageGeneric.X_NEG;
            	else if (view.equals("z_pos"))
            		resliceView = ImageGeneric.Z_POS;
            	else if (view.equals("y_pos"))
            		resliceView = ImageGeneric.Y_POS;
            	else if (view.equals("x_pos"))
            		resliceView = ImageGeneric.X_POS;
            }
            
            if(cmdLine.hasOption(OBJECT_CMDS))
            	objectCommands = cmdLine.getOptionValues(OBJECT_CMDS);
            
            if (cmdLine.hasOption(SAMPLINGRATE)) {
                samplingRate = Float.parseFloat(cmdLine.getOptionValue(SAMPLINGRATE));
            }
            if (cmdLine.hasOption(CHIMERAPORT)) {
                chimeraPort = Integer.parseInt(cmdLine.getOptionValue(CHIMERAPORT));
            }
            inverty = cmdLine.hasOption(INVERTY);
           
           
           
        } catch (Exception ex) {
        	ex.printStackTrace();
        }
    }
    
    public String getRenderLabel()
    {
        return renderLabels[0];
    }
    
    public String[] getDisplayLabels()
    {
        return displayLabels;
    }
    
    public void setBlock(String block)
    {
        this.block = block;
    }
    
    public String getBlock()
    {
        return block;
    }

    public boolean isMask() {
    	if(cmdLine == null)
    		return false;
        return cmdLine.hasOption("mask");
    }
   
    public Float getSamplingRate()
    {
        return samplingRate;
    }
    
    public Integer getChimeraPort()
    {
        return chimeraPort;
    }

    public void setChimeraPort(int i) {
        chimeraPort = i;
    }
   
}
