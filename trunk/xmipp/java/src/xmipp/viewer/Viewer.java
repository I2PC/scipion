package xmipp.viewer;

import ij.ImagePlus;
import xmipp.ij.XmippImageConverter;


public class Viewer {
	public static void main(String args[]) {
		Parameters parameters = new Parameters(args);
		try{
		for (String filename: parameters.files) {
			ImagePlus imp = XmippImageConverter.loadImage(filename);
			imp.getWindow();
			imp.show();
		}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
}



/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
//public class XmippBrowser implements PlugIn {
//
//    // Browser
//    String DIR;
//    boolean SINGLE_SELECTION = false;
//    String FILTER = "";
//    // Show
//    String INPUT[];
//    String MODE = COMMAND_PARAMETERS.MODE_DEFAULT;
//    boolean POLL = false;
//    boolean RENDER_IMAGES = false;
//    final static int MODE_TYPE_DEFAULT = 0;
//    final static int MODE_TYPE_IMAGE = 1;
//    final static int MODE_TYPE_TABLE = 2;
//    final static int MODE_TYPE_METADATA = 3;
//    int ROWS = -1, COLUMNS = -1;
//    int ZOOM = 1;

//    @Override
//    public void run(String string) {
//        Parameters parameters;
//
//        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
//            // "string" is used when called from another plugin or installed command.
//            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
//            parameters = new Parameters(Macro.getOptions().split(" "));
////            processArgs(Macro.getOptions().trim());
//        } else {    // From menu.
//            parameters = new Parameters();
//            parameters.directory = System.getProperty("user.dir");
//        }
//
//        if (parameters.files != null) {
//            openFiles(parameters);
//        }
//
//        if (parameters.directory != null) {
//            runBrowser(parameters);
//        }
//    }
//
//    void runBrowser(Parameters parameters) {
//        JDialogXmippBrowser frameBrowser = new JDialogXmippBrowser(
//                parameters.directory, parameters);
//        frameBrowser.setVisible(true);
//    }
//
//    static void openFiles(Parameters parameters) {
//    	String mode = parameters.mode.trim().toLowerCase();
//        if (mode.compareTo(COMMAND_PARAMETERS.OPENING_MODE_DEFAULT) == 0) {
//            ImagesWindowFactory.openFilesAsDefault(parameters.files, parameters);
//        } else if (mode.compareTo(COMMAND_PARAMETERS.OPENING_MODE_IMAGE) == 0) {
//            ImagesWindowFactory.openFilesAsImages(parameters.files, parameters);
//        } else if (mode.compareTo(COMMAND_PARAMETERS.OPENING_MODE_GALLERY) == 0) {
//            // Separates single images (index 0) from the rest of files (index 1).
//            String filesList[][] = getSeparatedFiles(parameters.files);
//
//            // Single images are opened in the same table...
//            if (filesList[0] != null) {
//                ImagesWindowFactory.openFilesAsGallery(filesList[0], true, parameters);
//            }
//
//            // ...while rest of files are opened in separated ones.
//            if (filesList[1] != null) {
//                ImagesWindowFactory.openFilesAsGallery(filesList[1], true, parameters);
//            }
//        } else if (parameters.mode.trim().toLowerCase().compareTo(COMMAND_PARAMETERS.OPENING_MODE_METADATA) == 0) {
//            ImagesWindowFactory.openFilesAsMetadata(parameters.files, parameters);
//        }
//    }
//
//    static String[][] getSeparatedFiles(String items[]) {
//        LinkedList<String> images = new LinkedList<String>();
//        LinkedList<String> files = new LinkedList<String>();
//        LinkedList<String> storeTo;
//
//        for (int i = 0; i < items.length; i++) {
//            String filename = items[i];
//
//            try {
//                storeTo = Filename.isSingleImage(filename) ? images : files;
//                storeTo.add(filename);
//            } catch (Exception e) {
//                IJ.error(e.getMessage());
//            }
//        }
//
//        String result[][] = new String[2][];
//
//        if (!images.isEmpty()) {
//            result[0] = images.toArray(new String[images.size()]);
//        }
//
//        if (!files.isEmpty()) {
//            result[1] = files.toArray(new String[files.size()]);
//        }
//
//        return result;
//    }
//
//    void processArgs(String args) {
//        String argsList[] = args.split(" ");
//        Options options = new Options();
//
//        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR, true, "");
//        options.addOption(COMMAND_PARAMETERS.OPTION_FILTER, true, "");
//        options.addOption(COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION, false, "");
//
//        options.addOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE, true, "");
//        options.addOption(COMMAND_PARAMETERS.OPTION_MODE, true, "");
//        options.addOption(COMMAND_PARAMETERS.OPTION_POLL, false, "");
//        options.addOption(COMMAND_PARAMETERS.OPTION_ZOOM, true, "");
//        options.addOption(COMMAND_PARAMETERS.OPTION_RENDER_IMAGES, false, "");
//        options.addOption(COMMAND_PARAMETERS.OPTION_TABLE_ROWS, true, "");
//        options.addOption(COMMAND_PARAMETERS.OPTION_TABLE_COLUMNS, true, "");
//
//        // It should be able to handle multiple files.
//        options.getOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE).setOptionalArg(true);
//        options.getOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE).setArgs(Integer.MAX_VALUE);
//
////        options.getOption(COMMAND_PARAMETERS.OPTION_MODE).setOptionalArg(true);
////        options.getOption(COMMAND_PARAMETERS.OPTION_MODE).setArgs(Integer.MAX_VALUE);
//
//        try {
//            BasicParser parser = new BasicParser();
//            CommandLine cmdLine = parser.parse(options, argsList);
//
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_DIR)) {
//                DIR = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_INPUT_DIR);
//            }
//
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_FILTER)) {
//                String filters[] = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_FILTER).split(COMMAND_PARAMETERS.FILTERS_SEPARATOR);
//
//                FILTER = "";
//                for (int i = 0; i < filters.length; i++) {
//                    FILTER += filters[i];
//                    if (i < filters.length - 1) {
//                        FILTER += " ";
//                    }
//                }
//            }
//
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_SINGLE_SELECTION)) {
//                SINGLE_SELECTION = true;
//            }
//
//            // Input.
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_INPUT_FILE)) {
//                INPUT = cmdLine.getOptionValues(COMMAND_PARAMETERS.OPTION_INPUT_FILE);
//            }
//
//            // Mode.
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_MODE)) {
//                MODE = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_MODE);
//            }
//
//            // Poll.
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_POLL)) {
//                POLL = true;
//            }
//
//            // Zoom.
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_POLL)) {
//                String str = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_ZOOM);
//                ZOOM = Integer.valueOf(str).intValue();
//            }
//
//            // Input.
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_RENDER_IMAGES)) {
//                RENDER_IMAGES = true;
//            }
//
//            // Table height.
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_TABLE_ROWS)) {
//                String str = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_TABLE_ROWS);
//                ROWS = Integer.valueOf(str).intValue();
//            }
//
//            // Table width.
//            if (cmdLine.hasOption(COMMAND_PARAMETERS.OPTION_TABLE_COLUMNS)) {
//                String str = cmdLine.getOptionValue(COMMAND_PARAMETERS.OPTION_TABLE_COLUMNS);
//                COLUMNS = Integer.valueOf(str).intValue();
//            }
//        } catch (Exception ex) {
//            throw new RuntimeException(ex);
//        }
//    }
//}
