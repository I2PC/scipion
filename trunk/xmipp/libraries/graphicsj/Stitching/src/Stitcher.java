
import ij.IJ;
import ij.ImagePlus;
import ij.io.Opener;
import ini.trakem2.ControlWindow;
import ini.trakem2.Project;
import ini.trakem2.display.Layer;
import ini.trakem2.display.LayerSet;
import ini.trakem2.display.Patch;
import java.awt.Rectangle;
import java.io.File;
import java.util.LinkedList;
import java.util.List;
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

    String filenames[];
    int x, y;
    String propertiesFile;
    String outputfilename;
    Parameters parameters;

    public Stitcher(String filenames[], int x, int y, String propertiesFile, String outputfilename) {
        this.filenames = fixFilenamesPaths(filenames);
        this.x = x;
        this.y = y;
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

//        double x = 0; // initial coordinates
//        double y = 0;

        for (int i = 0; i < filenames.length; i++) {
            ImagePlus imp = opener.openImage(filenames[i]);
            Patch patch = project.getLoader().addNewImage(imp, x, y);
            layer.add(patch);
            patches.add(patch);
        }

        // 5 - Register/Montage images
        //List<Patch> fixedPatches = new LinkedList<Patch>();

        alignPatches(patches, parameters);

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
/*        double scale = 1.0;
        int type = ImagePlus.GRAY8; // or .COLOR_RGB
        int c_alphas = 0xffffffff; // channel opacities, applies to RGB images only
        boolean quality = true; // use max possible rendering quality
        ArrayList list = null; // null means all. You can define a subset of Displayable objects
        //  of that layer (such as Patch objects) to be rendered instead.
        ImagePlus snapshot = project.getLoader().getFlatImage(layer, box, scale, c_alphas, type, Patch.class, list, quality);
        //new ij.io.FileSaver(snapshot).saveAsTiff("/Volumes/Data/jvega/Desktop/snapshot.tif");
        
        snapshot.setTitle("Montage [TrakEM2]");
        snapshot.show();*/

        // @TMP: Print results.
/*        System.out.println("----------------------");
        for (int i = 0; i < patches.size(); i++) {
        Patch patch = patches.get(i);
        Rectangle r = patch.getBoundingBox();
        AffineTransform t = patch.getAffineTransform();
        System.out.println(i + " > " + filenames[i] + ": ");
        System.out.println(" Position (x, y): (" + r.x + ", " + r.y + ")");
        System.out.println(" Affine transform: " + t);
        System.out.println("----------------------");
        }*/

//        System.err.println(" >>> Saving to: " + outputfilename);
//        ImagePlus ip = StackBuilder.buildStack(patches, box);
//        ip.setTitle("Montage [Stack]");
//        IJ.save(ip, outputfilename);
        if (StackBuilder.saveStack(patches, box, outputfilename)) {
            IJ.showMessage(" >>> File sucessfully saved: " + outputfilename);
        } else {
            System.err.println(" xxx ERROR: saving to: " + outputfilename);
        }

        ImagePlus imp = new ImagePlus(outputfilename);
        imp.show();
        /*        ImagePlus sum = StackBuilder.sumStack(ip);
        sum.setTitle("Montage [Sum]");
        sum.show();*/

        // 10 - Close the project
        project.destroy();
        // If GUI was disabled, be nice and reenable it:
        ControlWindow.setGUIEnabled(true);

        // Done!
        System.out.println("Done!!");
    }

    /**
     * Copied from @see mpicbg.trakem2.align.AlignTask to remove parameters window.
     * @param patches: the list of Patch instances to align, all belonging to the same Layer.
     * @param fixedPatches: the list of Patch instances to keep locked in place, if any.
     * @param parameters: parameters for alignment. If null, it will ask user for them.
     */
    public static void alignPatches(final List<Patch> patches, Parameters parameters) {
        AlignTask.alignPatches(parameters.getAlignParameters(),
                patches,
                null, //new LinkedList<Patch>(): fixed patches.
                true, // tilesAreInPlace
                false, // largestGraphOnly
                false, // hideDisconnectedTiles
                false // deleteDisconnectedTiles
                );
    }
}
