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
public class Param {

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
    public final static String OPENING_MODE_ROTSPECTRA = "rotspectra";
    public final static String RENDER_IMAGES = "render";
    public final static String SINGLE_SELECTION = "singlesel";
    public final static String SELECTION_TYPE = "seltype";
    public final static String SELECTION_TYPE_ANY = "any";
    public final static String SELECTION_TYPE_FILE = "file";
    public final static String SELECTION_TYPE_DIR = "dir";
    public final static String FILTER = "filter";
    public final static String FILTERS_SEPARATOR = ",";
    public final static String PORT = "port";
    public final static String DOWNSAMPLING = "downsampling";
    public final static String BAD_PIXELS_FACTOR = "factor";
    public final static String W1 = "w1";
    public final static String W2 = "w2";
    public final static String RAISED_W = "raised_w";
    public final static String MASKFILENAME = "mask";
    public final static String DEBUG = "debug";
    public final static String MASKTOOLBAR = "mask_toolbar";
    public final static String VIEW = "view";
    public final static String NO_GEO = "dont_apply_geo";
    public final static String NO_WRAP = "dont_wrap";
    
    
    
    public String directory;
    public String files[];
    public int port;
    public String filter;
    public boolean singleSelection;
    public String selectionType = SELECTION_TYPE_ANY;
    public String mode = OPENING_MODE_DEFAULT;
    public boolean poll;
    public int zoom = 0;
    public boolean renderImages = false;
    public String renderLabel = "first"; //Label to render, by default first
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

    
    
    

    public Param() {
    }

    public Param(String args[]) {
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
        options.addOption(SINGLE_SELECTION, false, "");
        options.addOption(SELECTION_TYPE, true, "");
        options.addOption(FILTER, true, "");

        options.addOption(OPENING_MODE, true, "");
        options.addOption(POLL, false, "");
        options.addOption(ZOOM, true, "");
        options.addOption(RENDER_IMAGES, true, "");
        options.addOption(DEBUG, false, "");
        options.addOption(MASKTOOLBAR, false, "");
        options.addOption(TABLE_ROWS, true, "");
        options.addOption(TABLE_COLUMNS, true, "");

        options.addOption(BAD_PIXELS_FACTOR, true, "");

        options.addOption(W1, true, "");
        options.addOption(W2, true, "");
        options.addOption(RAISED_W, true, "");

        options.addOption(DOWNSAMPLING, true, "");

        options.addOption(MASKFILENAME, true, "");

        options.addOption(INPUT_VECTORSFILE, true, "");
        options.addOption(INPUT_CLASSESFILE, true, "");
        options.addOption(INPUT_DATAFILE, true, "");
        options.addOption(VIEW, true, "z");
        options.addOption(NO_GEO, false, "");
        options.addOption(NO_WRAP, false, "");

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

            singleSelection = cmdLine.hasOption(SINGLE_SELECTION);

            if (cmdLine.hasOption(SELECTION_TYPE)) {
                selectionType = cmdLine.getOptionValue(SELECTION_TYPE);
            }

            if (cmdLine.hasOption(FILTER)) {
                String filters[] = cmdLine.getOptionValue(FILTER).split(FILTERS_SEPARATOR);

                filter = "";
                for (int i = 0; i < filters.length; i++) {
                    filter += filters[i];
                    if (i < filters.length - 1) {
                        filter += " ";
                    }
                }
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

            renderImages = cmdLine.hasOption(RENDER_IMAGES);
            if (renderImages){
            	renderLabel = cmdLine.getOptionValue(RENDER_IMAGES);
            }
            
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

            if (cmdLine.hasOption(MASKFILENAME)) {
                maskFilename = cmdLine.getOptionValue(MASKFILENAME);
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
            
           
        } catch (Exception ex) {
        	ex.printStackTrace();
        }
    }
}
