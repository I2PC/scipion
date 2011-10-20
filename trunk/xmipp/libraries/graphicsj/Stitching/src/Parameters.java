
import ij.IJ;
import java.io.FileInputStream;
import java.util.Properties;
import mpicbg.trakem2.align.Align;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Parameters {

    // Scale invariant interest point detector:
    final static String INITIAL_GAUSSIAN_BLUR = "initial_gaussian_blur";
    final static String STEPS_PER_SCALE_OCTAVE = "steps_per_scale_octave";
    final static String MINIMUM_IMAGE_SIZE = "minimum_image_size";
    final static String MAXIMUM_IMAGE_SIZE = "maximum_image_size";
    // Feature descriptor:
    final static String FEATURE_DESCRIPTOR_SIZE = "feature_descriptor_size";
    final static String FEATURE_DESCRIPTOR_ORIENTATION = "feature_descriptor_orientation_bins";
    final static String CLOSEST_NEXT_CLOSEST_RADIO = "closest/next_closest_ratio";
    // Geometric Consensus Filter:
    final static String MAXIMAL_ALIGNMENT_ERROR = "maximal_alignment_error";
    final static String MINIMUM_INLIER_RATIO = "minimum_inlier_ratio";
    final static String MINIMAL_NUMBER_OF_INLIERS = "minimal_number_of_inliers";
    final static String EXPECTED_TRANSFORMATION = "expected_transformation";// translation / rigid / similarity / affine
    final static String IGNORE_CONSTANT_BACKGROUND = "ignore_constant_background";
    final static String TOLERANCE = "tolerance";
    // Alignment:
    final static String DESIRED_TRANSFORMATION = "desire_transformation";
    final static String CORRESPONDENCE_WEIGTH = "correspondence_weight";
    final static String MAXIMAL_ITERATIONS = "maximal_iterations";
    final static String MAXIMAL_PLATEAU_WIDTH = "maximal_plateau_width";
    final static String FILTER_OUTLIERS = "filter_outliers";
    final static String MEAN_FACTOR = "mean_factor";
    final static String OVERLAP_MARGIN = "overlap_margin";
    // Miscellaneous:
    // tiles are roughly in place=yes
    // consider largest graph only=yes
    // hide tiles from non-largest graph=yes
    // delete tiles from non-largest graph=yes
    final static String[] modelStrings = new String[]{"Translation", "Rigid", "Similarity", "Affine"};
    Properties properties;

    /*    public static void main(String args[]) {
    String properties = "stitching.properties";
    Parameters parameters = new Parameters(properties);
    
    final Align.ParamOptimize p = parameters.getAlignParameters();
    final GenericDialog gd = new GenericDialog("Align Tiles");
    p.addFields(gd);
    
    gd.showDialog();
    }*/
    public Parameters(String filename) {
        properties = loadProperties(filename);
    }

    public Align.ParamOptimize getAlignParameters() {
        Align.ParamOptimize p = Align.paramOptimize.clone();

        // Sets parameters from properties.
        p.sift.initialSigma = getPropertyFloat(INITIAL_GAUSSIAN_BLUR);
        p.sift.steps = getPropertyInt(STEPS_PER_SCALE_OCTAVE);
        p.sift.minOctaveSize = getPropertyInt(MINIMUM_IMAGE_SIZE);
        p.sift.maxOctaveSize = getPropertyInt(MAXIMUM_IMAGE_SIZE);
        p.sift.fdSize = getPropertyInt(FEATURE_DESCRIPTOR_SIZE);
        p.sift.fdBins = getPropertyInt(FEATURE_DESCRIPTOR_ORIENTATION);

        p.rod = getPropertyFloat(CLOSEST_NEXT_CLOSEST_RADIO);
        p.maxEpsilon = getPropertyFloat(MAXIMAL_ALIGNMENT_ERROR);
        p.minInlierRatio = getPropertyFloat(MINIMUM_INLIER_RATIO);
        p.minNumInliers = getPropertyInt(MINIMAL_NUMBER_OF_INLIERS);
        p.expectedModelIndex = getModelIndex(modelStrings, getProperty(EXPECTED_TRANSFORMATION));
        p.rejectIdentity = getPropertyBoolean(IGNORE_CONSTANT_BACKGROUND);
        p.identityTolerance = getPropertyFloat(TOLERANCE);
        p.desiredModelIndex = getModelIndex(modelStrings, getProperty(DESIRED_TRANSFORMATION));
        p.correspondenceWeight = getPropertyFloat(CORRESPONDENCE_WEIGTH);

        p.filterOutliers = getPropertyBoolean(FILTER_OUTLIERS);
        p.meanFactor = getPropertyFloat(MEAN_FACTOR);

        return p;
    }

    public int getOverlapMargin() {
        return getPropertyInt(OVERLAP_MARGIN);
    }

    String getProperty(String property) {
        return properties.getProperty(property);
    }

    float getPropertyFloat(String property) {
        return Float.parseFloat(properties.getProperty(property));
    }

    double getPropertyDouble(String property) {
        return Double.parseDouble(properties.getProperty(property));
    }

    int getPropertyInt(String property) {
        return Integer.parseInt(properties.getProperty(property));
    }

    boolean getPropertyBoolean(String property) {
        return Boolean.parseBoolean(properties.getProperty(property));
    }

    int getModelIndex(String model[], String m) {
        for (int i = 0; i < model.length; i++) {
            if (m.toLowerCase().compareTo(model[i].toLowerCase()) == 0) {
                return i;
            }
        }

        return -1;
    }

    static Properties loadProperties(String filename) {
        Properties properties = new Properties();

        try {
            FileInputStream file = new FileInputStream(filename);

            properties.load(file);

            file.close();
        } catch (Exception ex) {
//            ex.printStackTrace();
            IJ.error(ex.getMessage());
        }

//        properties.list(System.out);
//        System.out.println("---------------------------------");

        return properties;
    }
}
