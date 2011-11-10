
import ij.ImagePlus;
import ij.io.Opener;
import ini.trakem2.ControlWindow;
import ini.trakem2.Project;
import ini.trakem2.display.Layer;
import ini.trakem2.display.LayerSet;
import ini.trakem2.display.Patch;
import java.awt.Rectangle;
import java.io.File;
import java.util.ArrayList;
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

    String filenames[];
    int xs[], ys[]; // Initial coordinates.
    String propertiesFile;
    String outputfilename;
    String stackfilename;
    Parameters parameters;

    public Stitcher(String filenames[], String propertiesFile, String outputfilename, int x, int y, String stackfilename) {
        this.filenames = fixFilenamesPaths(filenames);
        this.outputfilename = outputfilename;
        this.stackfilename = stackfilename;
        this.xs = new int[]{0, x};
        this.ys = new int[]{0, y};

        parameters = new Parameters(propertiesFile);
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
        LayerSet layerset = project.getRootLayerSet();
        layerset.add(new Layer(project, 1, 1, layerset));
        Layer layer = layerset.getLayer(0); // there is one layer by default in a new Project.  

        // 3 - Adding images to a Layer
        final Opener opener = new Opener();
        List<Patch> patches = new LinkedList<Patch>();

        for (int i = 0; i < filenames.length; i++) {
            ImagePlus imp = opener.openImage(filenames[i]);
            Patch patch = project.getLoader().addNewImage(imp, xs[i], ys[i]);
            layer.add(patch);
            patches.add(patch);
        }

        // 5 - Register/Montage images
        //List<Patch> fixedPatches = new LinkedList<Patch>();
        alignPatches(patches, parameters);

        // 9 - Take a snapshot of the registered images
        Rectangle box = layer.getMinimalBoundingBox(Patch.class, true); // or any other ROI you like
/*        double scale = 1.0;
        int type = ImagePlus.GRAY8; // or .COLOR_RGB
        int c_alphas = 0xffffffff; // channel opacities, applies to RGB images only
        boolean quality = true; // use max possible rendering quality
        ArrayList list = null; // null means all. You can define a subset of Displayable objects
        //  of that layer (such as Patch objects) to be rendered instead.
        ImagePlus snapshot = project.getLoader().getFlatImage(layer, box, scale, c_alphas, type, Patch.class, list, quality);
        new ij.io.FileSaver(snapshot).saveAsTiff("/home/juanjo/Desktop/snapshot.tif");

        snapshot.setTitle("Montage [TrakEM2]");
        snapshot.show();*/

        if (outputfilename != null) {
            if (ImagesBuilder.saveResult(patches, box, parameters.getOverlapMargin(), outputfilename)) {
                System.out.println(" >>> Results sucessfully saved: " + outputfilename);
            } else {
                System.err.println(" xxx ERROR: saving to: " + outputfilename);
            }
        }

        if (stackfilename != null) {
            if (ImagesBuilder.saveStack(patches, box, stackfilename)) {
                System.out.println(" >>> Stack sucessfully saved: " + stackfilename);
            } else {
                System.err.println(" xxx ERROR: saving to: " + stackfilename);
            }
        }

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
        // Fix first patch.
        List<Patch> fixed = new LinkedList<Patch>();
        fixed.add(patches.get(0));

        AlignTask.alignPatches(parameters.getAlignParameters(),
                patches,
                fixed, // fixed patches.
                true, // tilesAreInPlace
                false, // largestGraphOnly
                false, // hideDisconnectedTiles
                false // deleteDisconnectedTiles
                );
    }
}
