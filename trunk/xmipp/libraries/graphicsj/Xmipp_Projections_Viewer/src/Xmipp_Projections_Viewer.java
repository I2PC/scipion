
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
public class Xmipp_Projections_Viewer implements PlugIn, UniverseListener {

    static {
        // loads "XmippDataJava"
        System.loadLibrary("XmippDataJava");
    }
    private Sphere sphere;
    private Image3DUniverse universeVolume;
    private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
    private ImagePlus volumeIP, sphereIP;
    private MultidimArrayd xmippVolume;  // Volume for Xmipp library.
    private GlobalTransform gt = new GlobalTransform();
    private boolean dispatched = false;
    private ProjectionTimer timer = new ProjectionTimer(500);
    private static ProjectionWindow projectionWindow;
    private static JFrameImagesTable frameImagesTable;
    private final static String COMMAND_OPTION_VOLUME = "vol";
    private final static int INDEX_VOLUME = 0;

    public static void main(String args[]) {
        IJ.getInstance();
        (new Xmipp_Projections_Viewer()).run("");
    }

    public void run(String string) {
        String fileVolume = null;

        if (IJ.isMacro() && !Macro.getOptions().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            String argsList[] = processArgs(Macro.getOptions().trim());

            try {
                fileVolume = (new File(Macro.getOptions().trim())).getCanonicalPath();
            } catch (IOException ioex) {
                ioex.printStackTrace();
            }
            /*            if (!fileVolume.startsWith(File.separator)) {
            fileVolume = (new File(fileVolume)).getAbsolutePath();//System.getProperty("user.dir") + File.separator + fileVolume;
            }
            if (!fileEulerAngles.startsWith(File.separator)) {
            fileEulerAngles = System.getProperty("user.dir") + File.separator + fileEulerAngles;
            }*/
        } else {    // From menu.
            JFrameLoad frameLoad = new JFrameLoad();

            frameLoad.setLocationRelativeTo(null);
            frameLoad.setVisible(true);
            if (!frameLoad.cancelled) {
                fileVolume = frameLoad.getVolumeFile();
            }
        }
//        System.out.println("Volume: [" + fileVolume + "]");
//        System.out.println("Euler : [" + fileEulerAngles + "]");

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

            // Loads euler angles file...
            IJ.showStatus(LABELS.MESSAGE_LOADING_EULER_ANGLES_FILE);
            sphere = new Sphere();

            // ...and the sphere.
            IJ.showStatus(LABELS.MESSAGE_BUILDING_SPHERE);
            sphereIP = sphere.build();

            new StackConverter(sphereIP).convertToRGB();

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

    public static String[] processArgs(String args) {
        String argsList[] = args.split(" ");

        String parameters[] = {null, null};

        Options options = new Options();
        options.addOption(COMMAND_OPTION_VOLUME, true, "Volume file");

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Volume file.
            if (cmdLine.hasOption(COMMAND_OPTION_VOLUME)) {
                parameters[INDEX_VOLUME] = cmdLine.getOptionValue(COMMAND_OPTION_VOLUME);
                //System.out.println(" *** VOLUME: " + parameters[INDEX_VOLUME]);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return parameters;
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

        int projectionW = sphereIP.getWidth();
        int projectionH = sphereIP.getHeight();

        ImagePlus image = sphere.getProjection(xmippVolume, rot, tilt, projectionW, projectionH);
        int n_images = sphere.getFiles(rot, tilt).size();

        // Shows projection using a custom window.
        showProjection(image, n_images);
    }

    private void showProjection(final ImagePlus ip, final int n_images) {
        if (projectionWindow == null) {
            projectionWindow = new ProjectionWindow(ip);
        }

        projectionWindow.update(ip, n_images);
        projectionWindow.setVisible(true);
    }

    public void analyzeProjection() {
        double angles[] = getEyeAngles(universeVolume);

        int projectionW = sphereIP.getWidth();
        int projectionH = sphereIP.getHeight();

        String scoreFiles[] = sphere.analyzeProjection(xmippVolume, angles[0], angles[1], projectionW, projectionH);

        // Shows table with scores
        showScoreFiles(scoreFiles);

        // Let's try to clean up memory.
        System.gc();

        IJ.showStatus("");
    }

    private void showScoreFiles(String scoreFiles[]) {
        IJ.showStatus(LABELS.MESSAGE_LOADING_SCORE_FILE);

        if (frameImagesTable == null) {
            frameImagesTable = new JFrameImagesTable();
        }

        frameImagesTable.loadScoreFiles(scoreFiles);

        ImageWindow3D window = universeVolume.getWindow();
        Point location = window.getLocation();
        location.translate(0, window.getHeight());

        frameImagesTable.setLocation(location);

        frameImagesTable.setWidth(universeVolume.getWindow().getWidth());

        frameImagesTable.setVisible(true);
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
