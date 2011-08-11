
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ini.trakem2.ControlWindow;
import ini.trakem2.Project;
import ini.trakem2.display.Layer;
import ini.trakem2.display.LayerSet;
import ini.trakem2.display.Patch;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import mpicbg.trakem2.align.Align;
import mpicbg.trakem2.align.AlignTask;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Stitching implements PlugIn {

    private final static String OPTION_INPUT_FILE = "i";
    private final static String OPTION_INPUT_FILE_DESCRIPTION = "input file(s)";
    String FILES[];

    public static void main(String args[]) {
        String files[] = {
            "/Volumes/Data/jvega/Desktop/img-1.png",
            "/Volumes/Data/jvega/Desktop/img-2.png",
            "/Volumes/Data/jvega/Desktop/img-3.png",
            "/Volumes/Data/jvega/Desktop/img-4.png"
        };

        stitch(files);
    }

    public void run(String string) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            processArgs(Macro.getOptions().split(" "));

            stitch(FILES);
        } else {    // From menu.
            IJ.error("No menu implementation yet.");
        }
    }

    private void processArgs(String args[]) {
        Options options = new Options();

        options.addOption(OPTION_INPUT_FILE, true, OPTION_INPUT_FILE_DESCRIPTION);

        // It should be able to handle multiple files.
        options.getOption(OPTION_INPUT_FILE).setOptionalArg(true);
        options.getOption(OPTION_INPUT_FILE).setArgs(Integer.MAX_VALUE);

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, args);

            if (cmdLine.hasOption(OPTION_INPUT_FILE)) {
                FILES = cmdLine.getOptionValues(OPTION_INPUT_FILE);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }

    public static void stitch(String filenames[]) {
        // 0 - Disable GUI if headless desired
        ControlWindow.setGUIEnabled(false);

        // 1 - Create a new, empty project
        Project project = Project.newFSProject("blank", null, System.getProperty("java.io.tmpdir"));

        // 2 - Obtain LayerSet and Layer pointers
        LayerSet ls = project.getRootLayerSet();
        ls.add(new Layer(project, 1, 1, ls));
        Layer layer = ls.getLayer(0); // there is one layer by default in a new Project.  

        // 3 - Adding images to a Layer
        final Opener opener = new Opener();
        List<Patch> patches = new LinkedList<Patch>();

        double x = 0; // initial coordinates
        double y = 0;

        for (int i = 0; i < filenames.length; i++) {
            ImagePlus imp = opener.openImage(filenames[i]);
            Patch patch = project.getLoader().addNewImage(imp, x, y);
            layer.add(patch);
            patches.add(patch);
        }

        // 5 - Register/Montage images
        List<Patch> fixedPatches = new LinkedList<Patch>();

        alignPatches(patches, fixedPatches);

        // @TODO Check!
        // 6 - Adjust brightness and contrast

        // ... add all other images
        // WARNING: must all be of the same dimensions.
        //          Will refuse to work otherwise.
        //
        // Strategy A: set same min and max for all images,
        //   but shifting the histogram's peak to compensate a bit.
        double min = 15600; // for 16-bit TEM images
        double max = 24700;
        Thread task2 = project.getLoader().setMinAndMax(patches, min, max);
        if (task2 != null) {
            try {
                task2.join();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        // @TODO Try this approach!
        // Strategy B: order all tiles by the stdDev of their histograms, and then
        // use the central 50% to obtain average values to apply to all tiles.
/*        Thread task2 = project.getLoader().homogenizeContrast(al, null);
        if (null != task2) {
        try {
        task2.join();
        } catch (Exception e) {
        e.printStackTrace();
        }
        }*/

        // 9 - Take a snapshot of the registered images
        Rectangle box = layer.getMinimalBoundingBox(Patch.class); // or any other ROI you like
        double scale = 1.0;
        int type = ImagePlus.GRAY8; // or .COLOR_RGB
        int c_alphas = 0xffffffff; // channel opacities, applies to RGB images only
        boolean quality = true; // use max possible rendering quality
        ArrayList list = null; // null means all. You can define a subset of Displayable objects
        //  of that layer (such as Patch objects) to be rendered instead.
        ImagePlus snapshot = project.getLoader().getFlatImage(layer, box, scale, c_alphas, type, Patch.class, list, quality);
        //new ij.io.FileSaver(snapshot).saveAsTiff("/Volumes/Data/jvega/Desktop/snapshot.tif");

        // Done!
        System.out.println("Done!!");
        snapshot.show();

        // @TMP: Print results.
        System.out.println("----------------------");
        for (int i = 0; i < patches.size(); i++) {
            Patch patch = patches.get(i);
            Rectangle r = patch.getBoundingBox();
            AffineTransform t = patch.getAffineTransform();
            System.out.println(i + " > " + filenames[i] + ": ");
            System.out.println(" Position (x, y): (" + r.x + ", " + r.y + ")");
            System.out.println(" Affine transform: " + t);
            System.out.println("----------------------");
        }

        buildStack(patches, box);

        // 10 - Close the project
        project.destroy();
        // If GUI was disabled, be nice and reenable it:
        ControlWindow.setGUIEnabled(true);
    }

    private static void buildStack(List<Patch> patches, Rectangle box) {
        ImageStack is = new ImageStack(box.width, box.height);

        for (int i = 0; i < patches.size(); i++) {
            Patch patch = patches.get(i);

            try {
                float slice[] = buildSlice(patch, box);

                is.addSlice(patch.getImagePlus().getTitle(), slice);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        ImagePlus ip = new ImagePlus("Montage", is);
        ip.show();
    }

    private static float[] buildSlice(Patch patch, Rectangle box) throws Exception {
        ImagePlus ip = patch.getImagePlus();
        AffineTransform T = patch.getAffineTransform();
        T.invert(); // Inverts transformation.

        float pixels[] = (float[]) ip.getProcessor().convertToFloat().getPixels();
        Rectangle r = new Rectangle(ip.getWidth(), ip.getHeight());

        int w = box.width;
        int h = box.height;
        float slice[] = new float[w * h];

        Point p = new Point();

        for (int j = 0; j < h; j++) {
            p.y = j;
            for (int i = 0; i < w; i++) {
                p.x = i;

                Point p_ = new Point();
                T.transform(p, p_);

                // If point is inside image...
                if (r.contains(p_)) {

                    // Bilinear Interpolation.
                    float value = Interpolator.bilinear(pixels, ip.getWidth(), ip.getHeight(), p_.getX(), p_.getY());

                    // Store point.
                    slice[j * w + i] = value;//pixels[y * ip.getWidth() + x];
                }
            }
        }

        return slice;
    }

    // @TODO Parametrize.
    /**
     * Copied from @see mpicbg.trakem2.align.AlignTask to remove parameters window.
     * @param patches: the list of Patch instances to align, all belonging to the same Layer.
     * @param fixed: the list of Patch instances to keep locked in place, if any.
     */
    public static void alignPatches(final List< Patch> patches, final List< Patch> fixedPatches) {
        boolean tilesAreInPlace = false;
        boolean largestGraphOnly = false;
        boolean hideDisconnectedTiles = false;
        boolean deleteDisconnectedTiles = false;

//        if (patches.size() < 2) {
//            Utils.log("No images to align.");
//            return;
//        }
//
//        for (final Patch patch : fixedPatches) {
//            if (!patches.contains(patch)) {
//                Utils.log("The list of fixed patches contains at least one Patch not included in the list of patches to align!");
//                return;
//            }
//        }

        //final Align.ParamOptimize p = Align.paramOptimize;
//        final GenericDialog gd = new GenericDialog("Align Tiles");
//        Align.paramOptimize.addFields(gd);
//
//        gd.addMessage("Miscellaneous:");
//        gd.addCheckbox("tiles are rougly in place", tilesAreInPlace);
//        gd.addCheckbox("consider largest graph only", largestGraphOnly);
//        gd.addCheckbox("hide tiles from non-largest graph", hideDisconnectedTiles);
//        gd.addCheckbox("delete tiles from non-largest graph", deleteDisconnectedTiles);

//        gd.showDialog();
//        if (gd.wasCanceled()) {
//            return;
//        }

//        Align.paramOptimize.readFields(gd);
//        tilesAreInPlace = gd.getNextBoolean();
//        largestGraphOnly = gd.getNextBoolean();
//        hideDisconnectedTiles = gd.getNextBoolean();
//        deleteDisconnectedTiles = gd.getNextBoolean();

        Align.paramOptimize.sift.maxOctaveSize = 2024;
        final Align.ParamOptimize p = Align.paramOptimize.clone();

        AlignTask.alignPatches(p, patches, fixedPatches, tilesAreInPlace, largestGraphOnly, hideDisconnectedTiles, deleteDisconnectedTiles);
    }
}
