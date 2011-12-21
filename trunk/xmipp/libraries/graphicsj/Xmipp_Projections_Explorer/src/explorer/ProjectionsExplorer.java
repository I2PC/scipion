package explorer;

import sphere.Sphere;
import constants.LABELS;
import ij.IJ;
import ij.ImagePlus;
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
import java.text.DecimalFormat;
import java.util.ArrayList;
import javax.media.j3d.Transform3D;
import javax.media.j3d.View;
import javax.swing.Timer;
import javax.vecmath.Color3f;
import javax.vecmath.Point3d;
import sphere.Geometry;
import table.JFrameImagesTable;
import window.ProjectionWindow;
import xmipp.ImageGeneric;
import xmipp.Projection;
import xmippij.XmippImageConverter;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class ProjectionsExplorer implements UniverseListener {

    boolean use_sphere;
    Sphere sphere;
    Image3DUniverse universeVolume, universeSphere;
    final static int UNIVERSE_W = 400, UNIVERSE_H = 400;
    ImagePlus volumeIP, sphereIP;
    ImageGeneric xmippVolume;  // Volume for Xmipp library.
    GlobalTransform gt = new GlobalTransform();
    //boolean dispatched = false;
    ProjectionTimer timer = new ProjectionTimer(500);
    static ProjectionWindow projectionWindow;
    JFrameImagesTable frameImagesTable;
    final static String COMMAND_OPTION_INPUT = "i";
    final static String COMMAND_OPTION_EULER_ANGLES = "angles";
    final static int INDEX_VOLUME = 0;
    final static int INDEX_EULER_ANGLES = 1;

    public void run(String volumeFile, String eulerAnglesFile) {
        // Volume: required, Angles: optional.
        if (volumeFile != null && fileExists(volumeFile)) {
            if (eulerAnglesFile != null) {
                use_sphere = true;
                if (!fileExists(volumeFile)) {
                    IJ.error(LABELS.MESSAGE_ERROR_FNF(volumeFile));
                    return;
                }
            }
        } else {
            IJ.error(LABELS.MESSAGE_ERROR_FNF(volumeFile));
            return;
        }

        try {
            // Loads image plus for volume.
            IJ.showStatus(LABELS.MESSAGE_LOADING_VOLUME);
            xmippVolume = new ImageGeneric(volumeFile);
            xmippVolume.read(ImageGeneric.ALL_IMAGES);
            xmippVolume.setXmippOrigin();

            // Converts it to show.
            IJ.showStatus(LABELS.MESSAGE_CONVERTING_XMIPP_VOLUME);
            volumeIP = XmippImageConverter.convertImageGenericToImageJ(xmippVolume);
            //volumeIP.setTitle(volumeFile);
            //volumeIP.show();

            // Loads euler angles file...
            IJ.showStatus(LABELS.MESSAGE_LOADING_EULER_ANGLES_FILE);
            sphere = use_sphere ? new Sphere(eulerAnglesFile) : new Sphere();

            // ...and the sphere.
            if (use_sphere) {
                IJ.showStatus(LABELS.MESSAGE_BUILDING_SPHERE);
                sphereIP = sphere.build(xmippVolume);
                //sphereIP.show();
            }

            // Creates both universes and shows them.
            IJ.showStatus(LABELS.MESSAGE_BUILDING_VOLUME_UNIVERSE);
            int threshold = (int) getThresholdValue(xmippVolume);
            universeVolume = createVolumeUniverse(volumeIP, threshold);

            ImageWindow3D windowVolume = universeVolume.getWindow();
            windowVolume.setTitle(LABELS.TITLE_VOLUME);

            // Stores volume window location to place windows later.
            Point loc = windowVolume.getLocation();
            loc.translate(windowVolume.getWidth(), 0);

            if (use_sphere) {
                IJ.showStatus(LABELS.MESSAGE_BUILDING_SPHERE_UNIVERSE);
                universeSphere = createSphereUniverse(sphereIP);

                ImageWindow3D windowSphere = universeSphere.getWindow();

                windowSphere.setTitle(LABELS.TITLE_SPHERE);

                // Places windows one next to the other.
                windowSphere.setLocation(loc);

                // Sets location for projections window.
                loc.translate(windowSphere.getWidth(), 0);
            }

            // Starts by retrieving first projection.
            showProjection();

            // Places projection window next to sphere universe (or volume if the sphere is not being used).
            projectionWindow.setLocation(loc);

            IJ.showStatus(LABELS.MESSAGE_DONE);
        } catch (Exception ex) {
            IJ.error("Error at run: " + ex.getMessage());
            ex.printStackTrace();
        }
    }

    static double getThresholdValue(ImageGeneric image) {
        double threshold;
        try {
            double stats[] = image.getStatistics();
            double m = stats[0];
            double M = stats[1];
            double oldRange = M - m;

            double otsu = Projection.entropyOtsuSegmentation(image, 0.005, false);

            // Scales range.
            threshold = (otsu - m) * (255 / oldRange);    // value = ((oldValue - oldMin) * newRange / oldRange) + newMin
        } catch (Exception ex) {
            threshold = 255;
            IJ.error("ERROR getting otsu threshold value. Using: " + threshold);
        }

        return threshold;
    }

    private Image3DUniverse createSphereUniverse(ImagePlus imp) {
        Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W, UNIVERSE_H);
        universe.show();    // Shows...

        // Adds the image plus to universe.
        new StackConverter(imp).convertToRGB();
        Content c = universe.addContent(imp, Content.VOLUME);

        universe.select(c);
        c.setLocked(true);  // To avoid selected content independent rotations.

        // Adds listener to synchronize both universes.
        universe.addUniverseListener(this);

        return universe;
    }

    private Image3DUniverse createVolumeUniverse(ImagePlus imp, int threshold) {
        Image3DUniverse universe = new Image3DUniverse(UNIVERSE_W, UNIVERSE_H);

        // Adds the image plus to universe.
        new StackConverter(imp).convertToRGB();
        Content c = universe.addSurfacePlot(imp);/*), new Color3f(1f, 165f / 255, 82f / 255), imp.getTitle(),
        50, new boolean[]{true, true, true}, 1);*/
        c.displayAs(Content.SURFACE);
        c.setColor(new Color3f(1f, 165f / 255, 82f / 255));
        c.setName("test");
        universe.select(c);
        c.setLocked(true);  // To avoid selected content independent rotations.

        // Adds listener to synchronize both universes.
        universe.addUniverseListener(this);

        universe.show();    // Shows...

        return universe;
    }

    private synchronized void showProjection() throws Exception {
        double angles[] = getEyeAngles(universeVolume);

        double rot = angles[0];
        double tilt = angles[1];
        double psi = 0;//angles[2];

        ImageGeneric projection = sphere.getProjection(xmippVolume, getMatrix(universeVolume));
        ImagePlus imp = XmippImageConverter.convertImageGenericToImageJ(projection);
        imp.setTitle("Projection");

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
        try {
            double angles[] = getEyeAngles(universeVolume);

            int projectionW = sphereIP.getWidth();
            int projectionH = sphereIP.getHeight();

            String scoreFile = sphere.analyzeProjection(xmippVolume, getMatrix(universeVolume), projectionW, projectionH);

            // Shows table with scores
            showScoreFiles(scoreFile);
        } catch (Exception ex) {
            ex.printStackTrace();
            IJ.error(ex.getMessage());
        }

        // Let's try to clean up memory.
        System.gc();

        IJ.showStatus("");
    }

    private void showScoreFiles(String scoreFile) {
        if (!fileExists(scoreFile)) {
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

        return f.exists();
    }

    // @TODO Kino: Returns transformation matrix to retrieve projection.
    private double[] getMatrix(Image3DUniverse universe) {
        DecimalFormat df = new DecimalFormat("#.##");

        Transform3D t3d = new Transform3D();
        t3d.set(new double[]{
                    1, 0, 0, 0,
                    0, -1, 0, 0,
                    0, 0, -1, 0,
                    0, 0, 0, 1,});
        //t3d.rotX(Math.toRadians(180));

        Transform3D t = new Transform3D();
        universe.getRotationTG().getTransform(t);

        t3d.mul(t);

//        universe.getContent("").setTransform(new double[]{
//                    1, 0, 0, 0,
//                    0, -1, 0, 128,
//                    0, 0, -1, 128,
//                    0, 0, 0, 1,});

        double matrix[] = new double[4 * 4];
        t3d.get(matrix);

        System.out.println("M:");
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                double d = matrix[i * 4 + j];
                System.out.print(df.format(d) + " ");
            }
            System.out.println();
        }
        /*
        t3d.invert();
        t3d.get(matrix);
        //matrix[5] *= -1;
        matrix[10] *= -1;
        
        System.out.println("I:");
        for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
        double d = matrix[i * 4 + j];
        System.out.print(df.format(d) + " ");
        }
        System.out.println();
        }
        //System.out.println(matrix[2] + " / " + -matrix[6] + " / " + -matrix[10]);*/

        return matrix;
    }

    /**
     * Returns rotation and tilt for the current point of view.
     * @return
     */
    private double[] getEyeAngles(Image3DUniverse universe) throws Exception {
        Transform3D t3d = new Transform3D();

        universe.getCanvas().getView().getUserHeadToVworld(t3d);

        double matrix[] = new double[4 * 4];
        t3d.get(matrix);

        Point3d point = new Point3d(matrix[2], -matrix[6], -matrix[10]);

        //return Projection.eulerMatrix2Angles(matrix);
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

//        dispatched = false;
    }

    public void transformationStarted(View view) {
//        if (!dispatched) {  // Avoids repeated actions.
        timer.stop();
//        }
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

        //if (!dispatched) {  // Avoids repeated actions.
        timer.start();
        //}
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
            try {
                showProjection();
            } catch (Exception ex) {
                IJ.error("ERROR retrieving projection: " + ex.getMessage());
            }
        }
    }
}
