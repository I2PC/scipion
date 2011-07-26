
import ij.ImagePlus;
import ij.Macro;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ini.trakem2.ControlWindow;
import ini.trakem2.Project;
import ini.trakem2.display.Layer;
import ini.trakem2.display.LayerSet;
import ini.trakem2.display.Patch;
import java.awt.Rectangle;
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
public class Stitching implements PlugIn {

    public static void main(String args[]) {
        //new ImageJ();
        String filename1 = args[0];//"/Volumes/Data/jvega/Desktop/l.tif";
        String filename2 = args[1];//"/Volumes/Data/jvega/Desktop/r.tif";

        Stitching s = new Stitching();
        s.run(filename1 + " " + filename2);
    }

    public void run(String string) {
        String input[] = string.split(" ");

        String filename1 = input[0];
        String filename2 = input[1];

        stitch(filename1, filename2);
    }

    public static void stitch(String leftImg, String rightImg) {
        // Birth, life and death of a TrakEM2 project
        // ------------------------------------------
        // Create a project, add images to its default layer,
        // register and adjust the images,
        // export a snapshot of the result,
        // and save and close the project.

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
        ImagePlus imp1 = opener.openImage(leftImg);
        ImagePlus imp2 = opener.openImage(rightImg);
        double x = 0; // initial coordinates
        double y = 0;
        Patch patch1 = project.getLoader().addNewImage(imp1, x, y);
        layer.add(patch1);

        Patch patch2 = project.getLoader().addNewImage(imp2, x, y);
        layer.add(patch2);

        // Optionally, you can set the transform as well for each Patch:
        /*AffineTransform at = new AffineTransform();
        affine.translate(100, 300);
        affine.scale(0.5, 0.5);
        patch1.setAffineTransform(affine); // will copy its values,
        // not retain the AffineTransform object, so you can reuse it.*/

        // 4 - Update displays (optional, it's been done already at layer.add)
        //    You can use layer.addSilently(patch1) to avoid updating displays
        //    but such method is dangerous: no proper stack index ordering
        //    for different kinds of Displayable objects (such as text DLabel,
        //    which should be always on top, i.e. at the end of the layer list.)
        //Display.update(layer);

        // 5 - Register/Montage images
        List<Patch> patches = new LinkedList<Patch>();
        patches.add(patch1);
        patches.add(patch2);
        //for(int i=1;i<...;i++){
        //}
        // ... add all other images to the List
        // ...

        List<Patch> fixedPatches = new LinkedList<Patch>();

        Thread task = AlignTask.alignPatchesTask(patches, fixedPatches);
        // Automatic aligning: Hide parameters window.
        task.setName("Run$alignPatchesTask");
        Macro.setOptions(task, "max_size=2048");

        if (task != null) {
            try {
                task.join();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        // 6 - Adjust brightness and contrast
        ArrayList<Patch> al = new ArrayList<Patch>();
        al.add(patch1);

        // ... add all other images
        // WARNING: must all be of the same dimensions.
        //          Will refuse to work otherwise.
        //
        // Strategy A: set same min and max for all images,
        //   but shifting the histogram's peak to compensate a bit.
        double min = 15600; // for 16-bit TEM images
        double max = 24700;
        Thread task2 = project.getLoader().setMinAndMax(al, min, max);
        if (task2 != null) {
            try {
                task2.join();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

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

        // 7 - Update displays if any:
//        Display.update(layer);

        // 8 - Save the project
        // It is recommended to store .xml inside the storage folder, so that
        // then image file paths are relative and thus the whole folder
        // can be moved between different computers trivially:
        //project.saveAs(project.getLoader().getStorageFolder() + "the_name.xml", true); // overwriting any existing xml file with that name

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

        // 10 - Close the project
        project.destroy();
        // If GUI was disabled, be nice and reenable it:
        ControlWindow.setGUIEnabled(true);

        // Done!
        System.out.println("Done!!");

        Rectangle r1 = patch1.getBoundingBox();
        Rectangle r2 = patch2.getBoundingBox();
        System.out.println("----------------------");
        System.out.println("Position 1: " + r1);
        System.out.println("Position 2: " + r2);
        System.out.println("Ax: " + (r2.x - r1.x));
        System.out.println("Ay: " + (r2.y - r1.y));
        System.out.println("----------------------");
        System.out.println("Transform1: " + patch1.getAffineTransform());
        System.out.println("Transform1: " + patch2.getAffineTransform());
        System.out.println("----------------------");

        snapshot.show();

        // Cleanup: remove reference to the Thread and its associated options  
        Macro.setOptions(task, null);
    }
//
//    final static public Bureaucrat alignPatchesTask(final List< Patch> patches, final List< Patch> fixedPatches) {
//        if (0 == patches.size()) {
//            Utils.log("Can't align zero patches.");
//            return null;
//        }
//        Worker worker = new Worker("Aligning images", false, true) {
//
//            public void run() {
//                startedWorking();
//                try {
//                    AlignTask.alignPatches(patches, fixedPatches);
//                    Display.repaint();
//                } catch (Throwable e) {
//                    IJError.print(e);
//                } finally {
//                    finishedWorking();
//                }
//            }
//
//            public void cleanup() {
//                patches.get(0).getLayer().getParent().undoOneStep();
//            }
//        };
//
//        return Bureaucrat.createAndStart(worker, patches.get(0).getProject());
//    }
}
/*

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Display display = Display.getFront();
// get a layer from which to grab patches (class Patch wraps an image)
Layer layer = display.getLayer();

// grab all patches from the layer
ArrayList<Patch> all = (ArrayList<Patch>) layer.getDisplayables(Patch.class);

// the "nail" image won't move
Patch nail = all.get(0);
ArrayList<Patch> nailed = new ArrayList<Patch>();
nailed.add(nail);

// start registration
Thread task = AlignTask.alignPatchesTask( all, nailed );
try { if (null != task) task.join(); } catch (Exception e) { e.printStackTrace(); }

// optional: resize display to fit all images in it
layer.getParent().setMinimumDimensions();
 */
