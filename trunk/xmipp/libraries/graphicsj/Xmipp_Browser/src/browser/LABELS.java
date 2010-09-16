/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser;

import browser.files.FileBrowser;
import browser.imageitems.ImageItem;

/**
 *
 * @author Juanjo Vega
 */
public class LABELS {

    public static final String APP_NAME = "XMipp Browser";
    public static final String TITLE_ABOUT = "About " + APP_NAME;
    public static final String MESAGE_ABOUT = "About " + APP_NAME;
    public final static String TITLE_PREVIEW = "Preview";
    public final static String TITLE_ERROR = "Memory error";
    public final static String TITLE_SEND2WINDOW = "Send image to window";
    public final static String TITLE_MAIN_WINDOW = "ImageJ/Xmipp Browser";
    public final static String TITLE_GO_TO_SLICE = "Go to slice";

    public final static String TITLE_TABLE_WINDOW_IMAGES(int items) {
        return "ImageJ/Xmipp: " + items + " images";
    }

    public final static String TITLE_TABLE_WINDOW_VOLUME(ImageItem item) {
        return item.getFile().getName() + ": " + item.nslices + " slices.";
    }

    public final static String TITLE_TABLE_WINDOW_SELFILE(ImageItem item) {
        return item.getFile().getName() + ": " + item.nslices + " images.";
    }
    /**
     * Buttons labels
     */
    public final static String BUTTON_PARENT_DIRECTORY = "Parent";
    public final static String BUTTON_REFRESH_DIRECTORY = "Refresh";
    public final static String BUTTON_SEND_IMAGE2WINDOW = "Send2Window";
    public final static String BUTTON_OK = "Ok";
    public final static String BUTTON_CANCEL = "Cancel";
    public final static String BUTTON_OPEN = "Open Selected";
    public final static String BUTTON_SEND2TABLE = "Send to table";
    public final static String BUTTON_SHOW_LABELS = "Show labels";
    public final static String BUTTON_SELECT_ALL = "Select All";
    public final static String BUTTON_INVERT_SELECTION = "Invert Selection";
    public final static String BUTTON_REFRESH = "Refresh";
    public final static String BUTTON_MEAN = "Mean";
    public final static String BUTTON_STD_DEVIATION = "Std. Deviation";
    public final static String BUTTON_NORMALIZE = "Normalize";
    /**
     * Labels for operations
     */
    public final static String OPERATION_SELECT_OPERATION = "-- Select Operation --";
    public final static String OPERATION_MENU_FLIP = "Flip";
    public final static String OPERATION_MENU_ROTATE = "Rotate";
    public final static String OPERATION_MENU_INFO = "Info";
    public final static String OPERATION_MENU_PROCESS = "Process";
    public final static String OPERATION_MENU_RESLICE = "Reslice";
    public final static String OPERATION_MENU_OPEN_AS = "Open As";
    public final static String OPERATION_FLIP_VERTICAL = "Flip Vertical";
    public final static String OPERATION_FLIP_HORIZONTAL = "Flip Horizontal";
    public final static String OPERATION_ROTATE_LEFT = "Rotate Left";
    public final static String OPERATION_ROTATE_RIGHT = "Rotate Right";
    public final static String OPERATION_HISTOGRAM = "Histogram";
    public final static String OPERATION_PLOT_PROFILE = "Plot Profile";
    public final static String OPERATION_MEASUREMENTS = "Measurements";
    public final static String OPERATION_BC = "Bright&Contrast";
    public final static String OPERATION_THRESHOLD = "Threshold";
    public final static String OPERATION_SUBTRACTBG = "Subtract Background";
    public final static String OPERATION_GAUSSIAN_BLUR = "Gaussian Blur";
    public final static String OPERATION_FFT = "FFT";
    public final static String OPERATION_FFT_BAND_PASS_FILTER = "FFT Band Pass Filter";
    public final static String OPERATION_RESLICE_TOP = "Reslice Top";
    public final static String OPERATION_RESLICE_RIGHT = "Reslice Right";
    public final static String OPERATION_OPEN_AS_3D = "3D";
    public final static String OPERATION_MENU_GO_TO = "Go to";
    public final static String OPERATION_MENU_PLUGINS = "Plugins";
    public final static String OPERATION_MENU_ANALYSIS = "Analysis";
    public final static String OPERATION_MENU_CONVERT = "Convert";
    public final static String OPERATION_MENU_THRESHOLD = "Threshold";
    public final static String OPERATION_MENU_OTHER = "Other";
    /**
     * Labels for browser
     */
    public final static String LABEL_FILTER = "Filter: ";
    public final static String LABEL_WIDTH = "Width: ";
    public final static String LABEL_HEIGHT = "Height: ";
    public final static String LABEL_NSLICES = "Slices: ";
    public final static String LABEL_LOAD_PREVIEW = "Load preview";
    public final static String LABEL_TABLE_DISABLE = "Disable";
    public final static String LABEL_TABLE_ENABLE = "Enable";
    public final static String LABEL_TABLE_ENABLE_ALL = "Enable All";
    public final static String LABEL_TABLE_DISABLE_ALL = "Disable All";
    public final static String LABEL_TABLE_SAVE_AS_IMAGES = "Save as images";
    public final static String LABEL_TABLE_SAVE_AS_STACK = "Save as stack";
    public final static String LABEL_TABLE_SAVE_AS_SELFILE = "Save as selfile";
    /**
     * Labels for images panel
     */
    public final static String LABEL_ZOOM = "Zoom (%): ";
    public final static String LABEL_ROWS = "Rows: ";
    public final static String LABEL_COLUMNS = "Columns: ";
    public final static String LABEL_GO2IMAGE = "Go to image: ";
    public final static String LABEL_FLIP = "Flip: ";
    public final static String LABEL_ROTATE = "Rotate: ";
    public final static String LABEL_INFO = "Info: ";
    public final static String LABEL_PROCESS = "Process: ";
    public final static String LABEL_RESLICE = "Reslice: ";
    public final static String LABEL_GO_TO_SLICE = "Go to slice: ";
    public final static String LABEL_OPEN_AS = "Open as: ";

    public final static String LABEL_FILES_SHOWING(int elements, int total) {
        return total != elements ? " Showing " + elements + " / " + total + " elements." : "";
    }

    public final static String MESSAGE_MEMORY_ERROR(long fileSize, long maxMemory) {
        return "File size (" + FileBrowser.getFileSizeString(fileSize) + ") is bigger than available memory (" + FileBrowser.getFileSizeString(maxMemory) + ")";
    }
    public final static String MESSAGE_NO_ITEMS_SELECTED = "No items selected (or not enabled)";
}
