package constants;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class LABELS {

    public final static String APP_TITLE = "Xmipp Projection Explorer";
    public final static String TITLE_VOLUME = LABELS.APP_TITLE + ": Volume";
    public final static String TITLE_SPHERE = LABELS.APP_TITLE + ": Sphere";
    public final static String TITLE_PROJECTION = "Projection";
    public final static String LABEL_VOLUME_FILE = "Volume file:";
    public final static String LABEL_ANGLES_FILE = "Angles file:";
    public final static String LABEL_AVERAGE_IMAGE = "Average";
    //public final static String LABEL_VOLUME = "Volume:";
    //public final static String LABEL_ANGLES = "Angles:";

    public static String LABEL_N_IMAGES(int nimages) {
        return nimages + " images.";
    }
    public final static String BUTTON_LOAD_FILE = "Load";
    public final static String BUTTON_OK = "Ok";
    public final static String BUTTON_CANCEL = "Cancel";
    public final static String BUTTON_ANALYZE = "Analyze";
    public final static String MESSAGE_LOADING_VOLUME = "Loading Volume...";
    public final static String MESSAGE_LOADING_EULER_ANGLES_FILE = "Loading euler angles file...";
    public final static String MESSAGE_CONVERTING_XMIPP_VOLUME = "Converting Xmipp Volume...";
    public final static String MESSAGE_BUILDING_SPHERE = "Building sphere...";
    public final static String MESSAGE_BUILDING_VOLUME_UNIVERSE = "Building volume universe...";
    public final static String MESSAGE_BUILDING_SPHERE_UNIVERSE = "Building sphere universe...";
    public final static String MESSAGE_RETRIEVING_PROJECTION = "Retrieving projection...";
    public final static String MESSAGE_LOADING_SCORE_FILE = "Loading score file...";
    public final static String MESSAGE_ANALYZING_PROJECTION = "Analyzing projection...";
    public final static String MESSAGE_DONE = "Done!";

    public static String MESSAGE_ERROR_FNF(String filename) {
        return "ERROR: File not found: " + filename;
    }
}
