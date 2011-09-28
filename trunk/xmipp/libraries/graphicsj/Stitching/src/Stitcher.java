
import ij.IJ;
//import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.Opener;
import ini.trakem2.ControlWindow;
import ini.trakem2.Project;
import ini.trakem2.display.Layer;
import ini.trakem2.display.LayerSet;
import ini.trakem2.display.Patch;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import mpicbg.trakem2.align.Align;
import mpicbg.trakem2.align.AlignTask;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Stitcher implements Runnable {
/*    static {
        System.setProperty("plugins.dir", "/gpfs/fs1/home/bioinfo/jvega/xmipp/external/imagej/plugins");
    }*/

    String propertiesFile;
    String filenames[];
    String outputfilename;
    Parameters parameters;

    public Stitcher(String filenames[], String propertiesFile, String outputfilename) {
        this.filenames = fixFilenamesPaths(filenames);
        this.outputfilename = outputfilename;

        if (propertiesFile != null) {
            parameters = new Parameters(propertiesFile);
        }
    }

    // This method is needed 
    static String[] fixFilenamesPaths(String filenames[]) {
        String fixed[] = new String[filenames.length];

        for (int i = 0; i < filenames.length; i++) {
            if (!filenames[i].startsWith(File.separator)) {
                fixed[i] = System.getProperty("user.dir") + File.separator;
            } else {
                fixed[i] = "";
            }

            fixed[i] += filenames[i];
        }

        return fixed;
    }

    public void run() {
/*	ImageJ ij = IJ.getInstance();
        if (ij == null) {
        	ij = new ImageJ();
        }*/

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
	    System.out.println(">>> Filename["+i+"]: "+filenames[i]);
            ImagePlus imp = opener.openImage(filenames[i]);
            Patch patch = project.getLoader().addNewImage(imp, x, y);
            layer.add(patch);
            patches.add(patch);
        }

        // 5 - Register/Montage images
        List<Patch> fixedPatches = new LinkedList<Patch>();

        alignPatches(patches, fixedPatches, parameters);

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
                IJ.log("Exception: " + e.getMessage());
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
//        ImagePlus snapshot = project.getLoader().getFlatImage(layer, box, scale, c_alphas, type, Patch.class, list, quality);
        //new ij.io.FileSaver(snapshot).saveAsTiff("/Volumes/Data/jvega/Desktop/snapshot.tif");

        // Done!
        System.out.println("Done!!");
//        snapshot.setTitle("Montage [TrakEM2]");
//        snapshot.show();

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

        ImagePlus ip = StackBuilder.buildStack(patches, box);
        ip.setTitle("Montage [Stack]");
        if (outputfilename != null) {
            IJ.save(ip, outputfilename);
        } else {
            ip.show();
        }

        /*        ImagePlus sum = StackBuilder.sumStack(ip);
        sum.setTitle("Montage [Sum]");
        sum.show();*/

        // 10 - Close the project
        project.destroy();
        // If GUI was disabled, be nice and reenable it:
        ControlWindow.setGUIEnabled(true);
    }

    /**
     * Copied from @see mpicbg.trakem2.align.AlignTask to remove parameters window.
     * @param patches: the list of Patch instances to align, all belonging to the same Layer.
     * @param fixedPatches: the list of Patch instances to keep locked in place, if any.
     * @param parameters: parameters for alignment. If null, it will ask user for them.
     */
    public static void alignPatches(final List<Patch> patches, final List<Patch> fixedPatches, Parameters parameters) {
        boolean tilesAreInPlace = false;
        boolean largestGraphOnly = false;
        boolean hideDisconnectedTiles = false;
        boolean deleteDisconnectedTiles = false;

        Align.ParamOptimize p;

        // If alignment parameters are provided, it runs in auto mode. Otherwise, it will ask user for them.
        if (parameters == null) {
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

            p = Align.paramOptimize.clone();
            GenericDialog gd = new GenericDialog("Align Tiles");
            p.addFields(gd);

            gd.addMessage("Miscellaneous:");
            gd.addCheckbox("tiles are rougly in place", tilesAreInPlace);
            gd.addCheckbox("consider largest graph only", largestGraphOnly);
            gd.addCheckbox("hide tiles from non-largest graph", hideDisconnectedTiles);
            gd.addCheckbox("delete tiles from non-largest graph", deleteDisconnectedTiles);

            gd.showDialog();
            if (gd.wasCanceled()) {
                return;
            }

            p.readFields(gd);
            tilesAreInPlace = gd.getNextBoolean();
            largestGraphOnly = gd.getNextBoolean();
            hideDisconnectedTiles = gd.getNextBoolean();
            deleteDisconnectedTiles = gd.getNextBoolean();
        } else {
            p = parameters.getAlignParameters();
        }

        AlignTask.alignPatches(p, patches, fixedPatches, tilesAreInPlace, largestGraphOnly, hideDisconnectedTiles, deleteDisconnectedTiles);
    }
}
