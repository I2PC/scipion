
import window.JFrameLoad;
import sphere.Geometry;
import sphere.Sphere;
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
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import javax.media.j3d.Transform3D;
import javax.media.j3d.View;
import javax.swing.Timer;
import javax.vecmath.Point3d;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import table.JFrameImagesTable;
import window.ProjectionWindow;
import window.iAnalyzer;
import xmipp.ImageGeneric_;
import xmipp.Projection;
import xmippij.XmippImageConverter_;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Xmipp_Projections_Explorer implements PlugIn, UniverseListener, iAnalyzer {

    private boolean use_sphere = false;
    private Sphere sphere;
    private Image3DUniverse universeVolume, universeSphere;
    private final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
    private ImagePlus volumeIP, sphereIP;
    private ImageGeneric_ xmippVolume;  // Volume for Xmipp library.
    private GlobalTransform gt = new GlobalTransform();
    private boolean dispatched = false;
    private ProjectionTimer timer = new ProjectionTimer(500);
    private static ProjectionWindow projectionWindow;
    private static JFrameImagesTable frameImagesTable;
    private final static String COMMAND_OPTION_INPUT = "i";
    private final static String COMMAND_OPTION_EULER_ANGLES = "angles";
    private final static int INDEX_VOLUME = 0;
    private final static int INDEX_EULER_ANGLES = 1;

    public static void main(String args[]) {
        IJ.getInstance();
        //(new Xmipp_Projections_Explorer()).run("");
        String fileVolume = "/home/jvega/Escritorio/imgs_Roberto/kk.vol";
        String fileEulerAngles = null;

        /*try {
        ImagePlus imp = XmippImageConverter.convertToImagej(fileVolume);
        imp.show();
        } catch (Exception ex) {
        ex.printStackTrace();
        }*/

        Xmipp_Projections_Explorer xpe = new Xmipp_Projections_Explorer();
        xpe.run(fileVolume, fileEulerAngles);
    }

    public void run(String string) {
        String fileVolume = null;
        String fileEulerAngles = null;
        boolean runGui = false;

        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            String argsList[] = processArgs(Macro.getOptions().trim());

            if (argsList[INDEX_VOLUME] != null) {
                try {
                    fileVolume = (new File(argsList[INDEX_VOLUME])).getAbsoluteFile().getCanonicalPath();
                } catch (IOException ioex) {
                    System.out.println("Error loading " + argsList[INDEX_VOLUME]);
                }
            } else {
                fileVolume = null;
            }

            if (argsList[INDEX_EULER_ANGLES] != null) {
                try {
                    fileEulerAngles = (new File(argsList[INDEX_EULER_ANGLES])).getAbsoluteFile().getCanonicalPath();
                } catch (IOException ioex) {
                    System.out.println("Error loading " + argsList[INDEX_EULER_ANGLES]);
                }
            } else {
                fileEulerAngles = null;
            }

            if (fileVolume == null) {
                runGui = true;
            } else {
                runGui = false;
            }
        } else {    // From menu.
            runGui = true;
        }

        if (runGui) {
            JFrameLoad frameLoad = new JFrameLoad();

            frameLoad.setLocationRelativeTo(null);
            frameLoad.setVisible(true);
            if (!frameLoad.cancelled) {
                fileVolume = frameLoad.getVolumeFile();
                fileEulerAngles = frameLoad.getEulerAnglesFile();
            }
        }

        // Volume: required, Angles: optional.
        if (fileEulerAngles != null) {
            use_sphere = true;
        }

        if (fileVolume != null) {
            run(fileVolume, fileEulerAngles);
        }
    }

    public static String[] processArgs(String args) {
        String argsList[] = args.split(" ");

        String parameters[] = {null, null};

        Options options = new Options();
        options.addOption(COMMAND_OPTION_INPUT, true, "Volume file");
        options.addOption(COMMAND_OPTION_EULER_ANGLES, true, "Euler angles file");

        try {
            BasicParser parser = new BasicParser();
            CommandLine cmdLine = parser.parse(options, argsList);

            // Volume file.
            if (cmdLine.hasOption(COMMAND_OPTION_INPUT)) {
                parameters[INDEX_VOLUME] = cmdLine.getOptionValue(COMMAND_OPTION_INPUT);
            }

            // Angles file.
            if (cmdLine.hasOption(COMMAND_OPTION_EULER_ANGLES)) {
                parameters[INDEX_EULER_ANGLES] = cmdLine.getOptionValue(COMMAND_OPTION_EULER_ANGLES);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }

        return parameters;
    }

    private static void checkFile(String file) {
        File f = new File(file);
        if (!f.exists()) {
            IJ.error("File [" + file + "] not found.");
        }
    }

    public void run(String volumeFile, String eulerAnglesFile) {
        checkFile(volumeFile);
        if (use_sphere) {
            checkFile(eulerAnglesFile);
        }

        try {
            // Loads image plus for volume.
            IJ.showStatus(LABELS.MESSAGE_LOADING_VOLUME);
            xmippVolume = new ImageGeneric_(volumeFile);
            xmippVolume.setXmippOrigin();

            // Converts it to show.
            IJ.showStatus(LABELS.MESSAGE_CONVERTING_XMIPP_VOLUME);
            volumeIP = XmippImageConverter_.convertToImagej(xmippVolume);
            //volumeIP.setTitle(volumeFile);

            new StackConverter(volumeIP).convertToRGB();

            // Loads euler angles file...
            IJ.showStatus(LABELS.MESSAGE_LOADING_EULER_ANGLES_FILE);
            sphere = use_sphere ? new Sphere(eulerAnglesFile) : new Sphere();

            // ...and the sphere.
            if (use_sphere) {
                IJ.showStatus(LABELS.MESSAGE_BUILDING_SPHERE);
                sphereIP = sphere.build();

                new StackConverter(sphereIP).convertToRGB();
            }

            // Creates both universes and shows them.
            IJ.showStatus(LABELS.MESSAGE_BUILDING_VOLUME_UNIVERSE);
            double threshold = Projection.entropyOtsuSegmentation(xmippVolume, 0.005, false);
            universeVolume = createUniverse(volumeIP, Content.SURFACE, threshold);

            ImageWindow3D windowVolume = universeVolume.getWindow();
            windowVolume.setTitle(LABELS.TITLE_VOLUME);

            // Stores volume window location to place windows later.
            Point loc = windowVolume.getLocation();
            loc.translate(windowVolume.getWidth(), 0);

            if (use_sphere) {
                IJ.showStatus(LABELS.MESSAGE_BUILDING_SPHERE_UNIVERSE);
                universeSphere = createUniverse(sphereIP, Content.VOLUME);

                ImageWindow3D windowSphere = universeSphere.getWindow();

                windowSphere.setTitle(LABELS.TITLE_SPHERE);

                // Places windows one next to the other.
                windowSphere.setLocation(loc);

                // Sets location for projections window.
                loc.translate(windowSphere.getWidth(), 0);
            }

            // Starts by retrieving first projection.
// TODO Uncomment            showProjection();

            // Places projection window next to sphere universe (or volume if the sphere is not being used).
// TODO Uncomment            projectionWindow.setLocation(loc);

            IJ.showStatus(LABELS.MESSAGE_DONE);
        } catch (Exception ex) {
            IJ.error("Error at run: " + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    private Image3DUniverse createUniverse(ImagePlus volume, int type) {
        return createUniverse(volume, type, 1);
    }

    private Image3DUniverse createUniverse(ImagePlus volume, int type, double threshold) {
        Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W, UNIVERSE_H);
        universe.show();    // Shows...

        // Adds the sphere image plus to universe.
        Content c = universe.addContent(volume, type);
        universe.select(c);

        c.setLocked(true);  // To avoid selected content independent rotations.
        System.out.println(" *** Setting THRESHOLD to: " + threshold);
        c.setThreshold((int) Math.round(threshold));

        // Adds listener to synchronize both universes.
        universe.addUniverseListener(this);

        return universe;
    }

    private void showProjection() {
        double angles[] = getEyeAngles(universeVolume);

        double rot = angles[0];
        double tilt = angles[1];

        Projection projection = sphere.getProjection(xmippVolume, rot, tilt);
        ImagePlus imp = XmippImageConverter_.convertToImagej(projection, "Projection");

        ArrayList<String> images = sphere.getFiles(rot, tilt);
        int n_images = images != null ? images.size() : 0;

        // Shows projection using a custom window.
        showProjection(imp, n_images);
    }

    private void showProjection(final ImagePlus ip, final int n_images) {
        if (projectionWindow == null) {
            projectionWindow = new ProjectionWindow(ip, n_images, this, use_sphere);
        }

        projectionWindow.update(ip, n_images);
        projectionWindow.setVisible(true);
    }

    public void analyzeProjection() {
        double angles[] = getEyeAngles(universeVolume);

        int projectionW = sphereIP.getWidth();
        int projectionH = sphereIP.getHeight();

        String scoreFile = sphere.analyzeProjection(xmippVolume, angles[0], angles[1], projectionW, projectionH);

        // Shows table with scores
        showScoreFiles(scoreFile);

        // Let's try to clean up memory.
        System.gc();

        IJ.showStatus("");
    }

    private void showScoreFiles(String scoreFile) {
        if (fileExists(scoreFile)) {
            IJ.error(scoreFile + ": Not found!!");
        } else {
            IJ.showStatus(LABELS.MESSAGE_LOADING_SCORE_FILE);

            if (frameImagesTable == null) {
                frameImagesTable = new JFrameImagesTable();
            }

            frameImagesTable.loadScoreFile(scoreFile);

            ImageWindow3D window = universeVolume.getWindow();
            Point location = window.getLocation();
            location.translate(0, window.getHeight());

            frameImagesTable.setLocation(location);

            frameImagesTable.setWidth(
                    universeVolume.getWindow().getWidth() + universeSphere.getWindow().getWidth());

            frameImagesTable.setVisible(true);
        }
    }

    private boolean fileExists(String filename) {
        File f = new File(filename);

        return filename == null || filename.trim().isEmpty() || !f.exists();
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

    private void synchronizeUniverses(View view) {
        if (use_sphere) {
            if (view == universeVolume.getCanvas().getView()) {
                universeVolume.getGlobalTransform(gt);
                universeSphere.setGlobalTransform(gt);
            } else {
                universeSphere.getGlobalTransform(gt);
                universeVolume.setGlobalTransform(gt);
            }
        }

        dispatched = false;
    }

    public void transformationStarted(View view) {
        if (!dispatched) {  // Avoids repeated actions.
            timer.stop();
        }
    }

    public void transformationUpdated(final View view) {
        // Synchronizes both universes in a different thread to avoid locks.
        new Thread(
                new Runnable() {

                    public void run() {
                        synchronizeUniverses(view);
                    }
                }).start();
    }

    public void transformationFinished(View view) {
        synchronizeUniverses(view);

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
// TODO Uncomment            showProjection();
        }
    }
}
