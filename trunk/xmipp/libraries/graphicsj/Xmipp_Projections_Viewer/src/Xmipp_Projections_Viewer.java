
import constants.LABELS;
import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.plugin.PlugIn;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.DefaultUniverse.GlobalTransform;
import ij3d.Image3DUniverse;
import ij3d.ImageWindow3D;
import ij3d.UniverseListener;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.media.j3d.Transform3D;
import javax.swing.Timer;
import javax.vecmath.Point3d;
import table.JFrameImagesTable;
import window.ProjectionWindow;
import java.io.File;
import java.io.IOException;
import javax.media.j3d.View;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Xmipp_Projections_Viewer implements PlugIn, UniverseListener {

    static {
        // loads "XmippDataJava"
        System.loadLibrary("XmippDataJava");
    }
    private Sphere sphere;
    private Image3DUniverse universeVolume;
    private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
    private ImagePlus volumeIP;
    private MultidimArrayd xmippVolume;  // Volume for Xmipp library.
    private GlobalTransform gt = new GlobalTransform();
    private boolean dispatched = false;
    private ProjectionTimer timer = new ProjectionTimer(500);
    private static ProjectionWindow projectionWindow;
    private static JFrameImagesTable frameImagesTable;

    public static void main(String args[]) {
        IJ.getInstance();
        (new Xmipp_Projections_Viewer()).run("");
    }

    public void run(String string) {
        String fileVolume = null;

        if (IJ.isMacro() && !(Macro.getOptions() == null || Macro.getOptions().isEmpty())) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            fileVolume = Macro.getOptions().trim();

            try {
                fileVolume = (new File(Macro.getOptions().trim())).getCanonicalPath();

                if (!fileVolume.startsWith(File.separator)) {
                    fileVolume = (new File(fileVolume)).getAbsolutePath();//System.getProperty("user.dir") + File.separator + fileVolume;
                }
            } catch (IOException ioex) {
                ioex.printStackTrace();
            }
        } else {    // From menu.
            JFrameLoadVolume frameLoad = new JFrameLoadVolume();

            frameLoad.setLocationRelativeTo(null);
            frameLoad.setVisible(true);
            if (!frameLoad.cancelled) {
                fileVolume = frameLoad.getVolumeFile();
            }
        }

        if (fileVolume != null) {
            run(fileVolume, "");
        }
    }

    public void run(String volumeFile, String null_) {
        try {
            // Loads image plus for volume.
            IJ.showStatus(LABELS.MESSAGE_LOADING_VOLUME);
            volumeIP = IJ.openImage(volumeFile);

            IJ.showStatus(LABELS.MESSAGE_BUILDING_XMIPP_VOLUME);
            xmippVolume = ImageConverter.imagej2xmipp(volumeIP);   // Builds Xmipp Volume just once.

            new StackConverter(volumeIP).convertToRGB();

            sphere = new Sphere();

            // Creates both universes and shows them.
            IJ.showStatus(LABELS.MESSAGE_BUILDING_VOLUME_UNIVERSE);
            universeVolume = createUniverse(volumeIP, Content.VOLUME);

            ImageWindow3D windowVolume = universeVolume.getWindow();

            windowVolume.setTitle(LABELS.TITLE_VOLUME);

            // Places windows one next to the other.
            Point loc = windowVolume.getLocation();
            loc.translate(windowVolume.getWidth(), 0);

            showProjection();    // Retrieves first projection.

            projectionWindow.setLocation(loc);
        } catch (Exception ex) {
            IJ.write(ex.getMessage());
            ex.printStackTrace();
        }
    }

    private Image3DUniverse createUniverse(ImagePlus volume, int type) {
        Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W, UNIVERSE_H);
        universe.show();    // Shows...

        // Adds the sphere image plus to universe.
        Content c = universe.addVoltex(volume);
        c.displayAs(type);

        // Adds listener to update projection.
        universe.addUniverseListener(this);

        return universe;
    }

    private void showProjection() {
        double angles[] = getEyeAngles(universeVolume);

        double rot = angles[0];
        double tilt = angles[1];

        int projectionW = volumeIP.getWidth();
        int projectionH = volumeIP.getHeight();

        ImagePlus image = sphere.getProjection(xmippVolume, rot, tilt, projectionW, projectionH);

        // Shows projection using a custom window.
        showProjection(image);
    }

    private void showProjection(final ImagePlus ip) {
        if (projectionWindow == null) {
            projectionWindow = new ProjectionWindow(ip);
        }

        projectionWindow.update(ip);
        projectionWindow.setVisible(true);
    }

    /**
     * Returns rotation and tilt for the current point of view.
     * @return
     */
    private double[] getEyeAngles(Image3DUniverse universe) {
        Transform3D t3d = new Transform3D();
        universe.getCanvas().getView().getUserHeadToVworld(t3d);

        double matrix[] = new double[4 * 4];
        t3d.get(matrix);

        Point3d point = new Point3d(matrix[2], -matrix[6], -matrix[10]);

        return Geometry.getAngles(point);
    }

    public void transformationStarted(View view) {
        if (!dispatched) {  // Avoids repeated actions.
            timer.stop();
        }
    }

    public void transformationUpdated(final View view) {
    }

    public void transformationFinished(View view) {
        if (!dispatched) {  // Avoids repeated actions.
            timer.start();
        }
    }

    public void contentAdded(Content c) {
    }

    public void contentRemoved(Content c) {
    }

    public void contentChanged(Content c) {
    }

    public void contentSelected(Content c) {
    }

    public void canvasResized() {
    }

    public void universeClosed() {
    }

    class ProjectionTimer implements ActionListener {

        private Timer timer;

        public ProjectionTimer(int delay) {
            timer = new Timer(delay, this);

            timer.setRepeats(false);    // Doesn't repeats.
        }

        public void start() {
            timer.start();
        }

        public void stop() {
            timer.stop();
        }

        public void actionPerformed(ActionEvent e) {
            showProjection();
        }
    }
}
