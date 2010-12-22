/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser;

import browser.files.FileBrowser;
import browser.imageitems.listitems.AbstractImageItem;

/**
 *
 * @author Juanjo Vega
 */
public class LABELS {

    public final static String APP_NAME = "XMipp Browser";
    public final static String TITLE_ABOUT = "About " + APP_NAME;
    public final static String MESAGE_ABOUT = "About " + APP_NAME;
    public final static String TITLE_PREVIEW = "Preview";
    public final static String TITLE_ERROR = "Memory error";
    public final static String TITLE_SEND2WINDOW = "Send image to window";
    public final static String TITLE_MAIN_WINDOW = "ImageJ/Xmipp Browser";
    public final static String TITLE_GO_TO_SLICE = "Go to slice";

    public final static String TITLE_TABLE_WINDOW_IMAGES(int items) {
        return "ImageJ/Xmipp: " + items + " images";
    }

    public final static String TITLE_TABLE_WINDOW_VOLUME(AbstractImageItem item) {
        return item.getFile().getName() + ": " + item.nslices + " slices.";
    }

    public final static String TITLE_TABLE_WINDOW_SELFILE(String file, int nslices) {
        return file + ": " + nslices + " images.";
    }
    /**
     * Buttons labels
     */
    public final static String BUTTON_TO_STACK = "To Stack";
    public final static String BUTTON_PARENT_DIRECTORY = "Parent";
    public final static String BUTTON_REFRESH_DIRECTORY = "Refresh";
    public final static String BUTTON_CAPTURE_WINDOW = "Capture";
    public final static String BUTTON_OK = "Ok";
    public final static String BUTTON_CANCEL = "Cancel";
    public final static String BUTTON_OPEN = "Open Selected";
    public final static String BUTTON_SEND2TABLE = "Send to table";
    public final static String BUTTON_SHOW_LABELS = "Show labels";
    public final static String BUTTON_SELECT_ALL = "Select All";
    public final static String BUTTON_INVERT_SELECTION = "Invert Selection";
    public final static String BUTTON_REFRESH = "Refresh";
    public final static String BUTTON_AVERAGE = "Average";
    public final static String BUTTON_STD_DEVIATION = "Std. Deviation";
    public final static String BUTTON_NORMALIZE = "Normalize";
    /**
     * Labels for operations
     */
    public final static String OPERATION_MENU_TRANSFORM = "Transform";
    public final static String OPERATION_TO_TABLE = "To Table";
    public final static String OPERATION_MENU_FLIP = "Flip";
    public final static String OPERATION_MENU_ROTATE = "Rotate";
    public final static String OPERATION_CROP = "Crop";
    public final static String OPERATION_MENU_DENOISING = "Denoising";
    public final static String OPERATION_MENU_OPEN_WITH = "Open with...";
    public final static String OPERATION_MENU_THRESHOLDING = "Thresholding";
    public final static String OPERATION_MENU_BINARY = "Binary";
    public final static String OPERATION_MENU_INFO = "Info";
    public final static String OPERATION_MENU_PROCESS = "Process";
    public final static String OPERATION_MENU_RESLICE = "Reslice";
    public final static String OPERATION_MENU_GO_TO = "Go to";
    public final static String OPERATION_MENU_THRESHOLD = "Threshold";
    public final static String OPERATION_MENU_DRAW = "Draw";
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
    public final static String OPERATION_STACK_SLICER = "StackSlicer (dynamic orthogonal views)";
    public final static String OPERATION_RADIAL_RESLICE = "Radial Reslice (orthogonal reconstructions)";
    public final static String OPERATION_STRAIGHTEN_CURVED_OBJECTS = "Straighten Curved Objects";
    public final static String OPERATION_UNTILT_STACK = "Untilt Stack (removes tilt from a stack of images)";
    public final static String OPERATION_TRANSFORMJ = "TransformJ";
    public final static String OPERATION_TJABOUT = "TJ About";
    public final static String OPERATION_TJAFFINE = "TJ Affine";
    public final static String OPERATION_TJCROP = "TJ Crop";
    public final static String OPERATION_TJEMBED = "TJ Embed";
    public final static String OPERATION_TJMIRROR = "TJ Mirror";
    public final static String OPERATION_TJOPTIONS = "TJ Options";
    public final static String OPERATION_TJPANEL = "TJ Panel";
    public final static String OPERATION_TJROTATE = "TJ Rotate";
    public final static String OPERATION_TJSCALE = "TJ Scale";
    public final static String OPERATION_TJSHIFT = "TJ Shift";
    public final static String OPERATION_TJTURN = "TJ Turn";
    public final static String OPERATION_TJWEBSITE = "TJ Website";
    public final static String OPERATION_SURFACEJ = "SurfaceJ";
    public final static String OPERATION_OPEN_AS_3D = "3D Viewer";
    public final static String OPERATION_DUPLICATE = "Duplicate";
    public final static String OPERATION_GO_TO_SLICE = "Go to slice";
    public final static String OPERATION_SUBSTACK_MAKER = "Substack Maker";
    public final static String OPERATION_ANISOTROPIC_DIFFUSION = "Anisotropic Diffusion";
    public final static String OPERATION_MEAN_SHIFT_FILTER = "Mean Shift Filter";
    public final static String OPERATION_VOLUME_VIEWER_3D_SLICER = "Volume Viewer/3D-Slicer";
    public final static String OPERATION_VOLUMEJ = "VolumeJ";
    public final static String OPERATION_AUTOCORRELATION = "Autocorrelation";
    public final static String OPERATION_CONCENTRIC_CIRCLES = "Concentric Circles";
    public final static String OPERATION_LINE_ANALYZER = "Line Analyzer";
    public final static String OPERATION_OVAL_PROFILER_PLOT = "Oval Profile Plot";
    public final static String OPERATION_RADIAL_PROFILE_PLOT_ANGLE = "Radial Profile Plot Angle";
    public final static String OPERATION_RADIAL_PROFILE_PLOT_HEIGHT = "Radial Profile Plot Height";
    public final static String OPERATION_CONTOUR_PLOTTER = "Contour Plotter";
    public final static String OPERATION_FEATUREJ = "FeatureJ";
    public final static String OPERATION_SYNC_WINDOWS = "Sync Windows";
    public final static String OPERATION_SYNC_MEASURE_3D = "Sync Measure 3D";
    public final static String OPERATION_OTSU_THRESHOLDING = "Otsu Thresholding";
    public final static String OPERATION_MULTI_OTSU_THRESHOLDING = "Multi Otsu Threshold";
    public final static String OPERATION_MIXTURE_MODELING_THRESHOLDING = "Mixture Modeling Thresholding";
    public final static String OPERATION_MAXIMUM_ENTROPY_THRESHOLDING = "Maximum Entropy Thresholding";
    public final static String OPERATION_MULTI_THRESHOLDER = "MultiThresholder";
    public final static String OPERATION_SIOX = "SIOX Simple Interactive Object Extraction";
    public final static String OPERATION_RATS = "RATS Robust Automatic Threshold Selection";
    public final static String OPERATION_VOXEL_COUNTER = "Voxel Counter";
    public final static String OPERATION_ERODE = "Erode";
    public final static String OPERATION_DILATE = "Dilate";
    public final static String OPERATION_OPEN = "Open";
    public final static String OPERATION_CLOSE = "Close";
    public final static String OPERATION_FLOAT_MORPHOLOGY = "Float Morphology";
    public final static String OPERATION_OUTLINE = "Outline";
    public final static String OPERATION_FILL_HOLES = "Fill holes";
    public final static String OPERATION_SKELETONIZE = "Skeletonize";
    public final static String OPERATION_DISTANCE_MAP = "Distance map";
    public final static String OPERATION_ULTIMATE_POINTS = "Ultimate points";
    public final static String OPERATION_WATERSHED = "Watershed";
    public final static String OPERATION_VORONOI = "Voronoi";
    public final static String OPERATION_OPTIONS = "Options";
    public final static String OPERATION_CLAHE = "CLAHE Contrast Limited Adaptive Histogram Equalization";
    public final static String OPERATION_PCA = "PCA Principal Component Analysis";
    public final static String OPERATION_POLAR_TRANSFORMER = "Polar Transformer convert to polar coordinates";
    public final static String OPERATION_DOTTED_AND_DASHED_LINES = "Dotted and Dashed Lines";
    public final static String OPERATION_RADIAL_GRID = "Radial Grid";
    public final static String OPERATION_TEMPLATE_MATCHING = "Template Matching";
    //public final static String OPERATION_FLOWJ = "FlowJ";
    //public final static String OPERATION_FLOW3J = "Flow3J";
    //public final static String OPERATION_OBJECTJ = "ObjectJ";
    /**
     * Labels for browser
     */
    public final static String LABEL_AUTO_AJUST_COLUMNS = "Auto adjust columns";
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
    public final static String LABEL_OPEN_AS = "Open as: ";

    public final static String LABEL_FILES_SHOWING(int elements, int total) {
        return total != elements ? " Showing " + elements + " / " + total + " elements." : "";
    }

    public final static String MESSAGE_MEMORY_ERROR(long fileSize, long maxMemory) {
        return "File size (" + FileBrowser.getFileSizeString(fileSize) + ") is bigger than available memory (" + FileBrowser.getFileSizeString(maxMemory) + ")";
    }
    public final static String MESSAGE_NO_ITEMS_SELECTED = "No items selected (or not enabled)";
}
