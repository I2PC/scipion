package sphere;

import constants.LABELS;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import xmipp.ImageGeneric;
import xmipp.Projection;
import xmipp.MDLabel;
import xmipp.MetaData;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Sphere {

    //private final static int R = 50;   // Sphere radius.
    private final static int ROT = Geometry.MAX_ROT;   // [0 - 359] = 360 elements.
    private final static int TILT = Geometry.MAX_TILT + 1;   // [0 - 180] = 181 elements.
    private final static int CRUST_THICKNESS = 2;
    private final static double MIN_SPREAD_VALUE = 1;   // Minimum value to spread after a hit.
    public final static int VALUE = 255;   // Initial value for hits.
    private SphereCell sphereMap[][];
    private final static String TEMPDIR_PATH = System.getProperty("java.io.tmpdir") + "/explorer";
    private File tempDir;
    private final static String XMIPP_CLASSIFY_ANALYZE_CLUSTER = "xmipp_classify_analyze_cluster";

    public Sphere() {
        sphereMap = new SphereCell[ROT][TILT];
        for (int i = 0; i < sphereMap.length; i++) {
            for (int j = 0; j < sphereMap[0].length; j++) {
                sphereMap[i][j] = new SphereCell();
            }
        }

    }

    public Sphere(String fileName) {
        this();
        loadAngles(fileName);

        equalizeValues(sphereMap);

        // Creates temporary work directory.
        tempDir = new File(TEMPDIR_PATH);
        try {
            Runtime.getRuntime().exec("rm -rf " + TEMPDIR_PATH);
            Runtime.getRuntime().exec("mkdir " + TEMPDIR_PATH);
        } catch (IOException ex) {
            IJ.error("Error creating sphere: " + ex.getMessage());
            throw new RuntimeException(ex);
        }

        tempDir.deleteOnExit();
    }

    private SphereCell getTableCell(Point2d point) {
        SphereCell cell = null;

        try {
            cell = sphereMap[(int) point.x][(int) point.y];
        } catch (Exception ex) {
        }

        return cell;
    }

    private void loadAngles(String filename) {
        try {
            IJ.showStatus(LABELS.MESSAGE_LOADING_EULER_ANGLES_FILE);
            MetaData md = new MetaData(filename);

            //File file = new File(filename);
            //String baseDir = file.getParent() + File.separator;

            long ids[] = md.findObjects();

            for (long id : ids) {
                String imageFileName = md.getValueString(MDLabel.MDL_IMAGE, id, true);

                //if (!imageFileName.startsWith("/")) {   // Builds absolute path (if it's not absolute already)
                //    imageFileName = baseDir + imageFileName;
                //}

                double rot = md.getValueDouble(MDLabel.MDL_ANGLEROT, id);
                double tilt = md.getValueDouble(MDLabel.MDL_ANGLETILT, id);

                double angles[] = Geometry.normalizeAngles(rot, tilt);
                Point2d hit = new Point2d(angles[0], angles[1]);

                spreadHit(hit, imageFileName);
            }
        } catch (Exception ex) {
            IJ.error("Error loading angles: " + ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    private void spreadHit(Point2d hitted, String imageFileName) {
        setLocked(hitted, true);   // Locks hitted point...

        double value = addHit(hitted, hitted, imageFileName);    // ...adds a hit...

        ArrayList<Point2d> points = new ArrayList<Point2d>();
        spreadHit(hitted, hitted, value, imageFileName, points);    // ...and spreads it...

        setLocked(hitted, false);   // Finally unlocks all the points.
        for (int i = 0; i < points.size(); i++) {
            setLocked(points.get(i), false);
        }
    }

    private void spreadHit(Point2d hitted, Point2d current, double value,
            String imageFileName, ArrayList<Point2d> points) {
        if (value > MIN_SPREAD_VALUE) {
            // Calculates all neighbours.
            ArrayList<Point2d> neighbours = getNeighbours(current);
            double values[] = new double[neighbours.size()];    // Individual values.

            // For each neighbour...
            for (int i = 0; i < neighbours.size(); i++) {
                Point2d neighbour = neighbours.get(i);

                // ...adds value (stores it to spread it later)...
                values[i] = addHit(neighbour, hitted, imageFileName);

                // ...and adds point to the locked list.
                points.add(neighbour);
            }

            // After every neighbour is locked an calculated, spreads it
            // recursively (keeping the last value added to know when to stop it.)
            for (int i = 0; i < neighbours.size(); i++) {
                spreadHit(hitted, neighbours.get(i), values[i], imageFileName, points);
            }
        }
    }

    /**
     * Gets all (not locked) neighbour points.
     * @param center
     * @return
     */
    private ArrayList<Point2d> getNeighbours(Point2d center) {
        ArrayList<Point2d> points = new ArrayList<Point2d>();

        // For all neighbours.
        for (double i = center.x - 1; i <= center.x + 1; i++) {
            for (double j = center.y - 1; j <= center.y + 1; j++) {
                Point2d point = new Point2d(i, j);

                Geometry.fixPoint(point);

                // If point is not yet locked...
                if (!isLocked(point)) {
                    points.add(point);  // ...stores and...
                    setLocked(point, true);    // ...locks it to avoid duplicates.
                }
            }
        }

        return points;
    }

    /**
     * Returns the new value added.
     * @param point
     * @param center
     * @param imageFileName
     * @return
     */
    private double addHit(Point2d point, Point2d center, String imageFileName) {
        SphereCell cell = getTableCell(point);

        // Adds a new fileName.
        cell.fileNames.add(imageFileName);

        // Gets value according to the 3D point related to rot and tilt (x and y at each point)
        double value = calculate3DPointValue(
                Geometry.getSphereCoordinates(center.x, center.y),
                Geometry.getSphereCoordinates(point.x, point.y), VALUE);

        cell.value += value;

        return value;
    }

    /**
     *
     * @return
     * value * e^(-1/2 * ((angle(p0,p1)/sigma)^2))
     */
    private static double calculate3DPointValue(Point3d p0, Point3d p1, double value) {
        double sigma = 0.05;

        // To vectors...
        Vector3d v0 = new Vector3d(p0);
        Vector3d v1 = new Vector3d(p1);

        // ...so getting its angle is easier.
        double angle = v0.angle(v1);

        return value * Math.pow(Math.E, -0.5 * Math.pow((angle / sigma), 2));
    }

    private void setLocked(Point2d point, boolean block) {
        getTableCell(point).locked = block;
    }

    private boolean isLocked(Point2d point) {
        return getTableCell(point).locked;
    }

    /**
     * Multiplies 
     * @param table
     */
    private static void equalizeValues(SphereCell table[][]) {
        // Gets minimum and maximum value in table...
        double minmax[] = getMinMaxValues(table);

        scaleValues(table, minmax[0], minmax[1]);   // ...and scales its values.
    }

    private static double[] getMinMaxValues(SphereCell table[][]) {
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;

        for (int i = 0; i < table.length; i++) {
            for (int j = 0; j < table[0].length; j++) {
                double value = table[i][j].value;
                if (value > max) {
                    max = value;
                }
                if (value < min) {
                    min = value;
                }
            }
        }

        return new double[]{min, max};
    }

    private static void scaleValues(SphereCell table[][], double min, double max) {
        double scale = VALUE / (max - min);

        for (int i = 0; i < table.length; i++) {
            for (int j = 0; j < table[0].length; j++) {
                //table[i][j].value *= scale;
                table[i][j].value = (table[i][j].value - min) * scale;
            }
        }
    }

    public ImagePlus build(ImageGeneric volume) throws Exception {
        int W, H, D;
        int radius = Math.min(Math.min(volume.getXDim() / 2, volume.getYDim() / 2), volume.getZDim() / 2);
        W = H = D = 2 * radius + 1;  // For each dimension: r - 0 - r = 2r + 1

        ImageProcessor ip = new ColorProcessor(W, H);
        ImagePlus sphere = new ImagePlus("Sphere", ip);

        // Room for output ~ [slice][image]
        int[][] stack = new int[D][W * H];

        buildSphere(stack, W, H, D, radius, CRUST_THICKNESS);

        // Builds the stack.
        ImageStack outputStack = new ImageStack(W, H, D);
        for (int slice = 0; slice < stack.length; slice++) {
            outputStack.setPixels(stack[slice], slice + 1);
        }

        sphere.setStack("Sphere", outputStack);

        return sphere;
    }

    private void buildSphere(int stack[][], int W, int H, int D, int r, int thickness) throws Exception {
        for (int z = -r; z < r; z++) {
            for (int x = -r; x < r; x++) {
                for (int y = -r; y < r; y++) {
                    Point3d point = new Point3d(x, -y, -z);   // (Translates point to center).

                    // Check point position related to sphere.
                    double distance = Geometry.length(point);//point.distance(center);

                    // Point is inside sphere.
                    if (distance <= r) {
                        int color = 0;

                        if (distance > r - thickness) {    // At sphere crust.
                            //double angles[] = Projection.eulerDirection2Angles(new double[]{point.x, point.y, point.z});
                            double angles[] = Geometry.getAngles(point);
                            angles = Geometry.absAngles(angles[0], angles[1]);

                            // Centers point.
                            point.x += W / D;
                            point.y += H / D;
                            point.z += D / D;

                            color = (int) Interpolator.bilinear(sphereMap, angles[0], angles[1]);   // ...to add the current value.
                        }

                        // Sets voxel color.
                        stack[z + D / 2][(y + H / 2) * W + (x + W / 2)] = ColorMap.getColor(color, VALUE + 1);
                    }
                }
            }
        }
    }

    public ImagePlus getTableAsImage() {
        int W = sphereMap.length;
        int H = sphereMap[0].length;

        // Output image.
        ImageProcessor ip = new ColorProcessor(W, H);
        ImagePlus image = new ImagePlus("Table", ip);

        int[] outPixels = (int[]) image.getProcessor().getPixels();

        for (int rot = 0; rot < sphereMap.length; rot++) {
            for (int tilt = 0; tilt < sphereMap[0].length; tilt++) {
                int color = (int) sphereMap[rot][tilt].value;  // Gets value.

                outPixels[tilt * W + rot] = ColorMap.getColor(color, VALUE + 1);
            }
        }

        // Sets pixels.
        image.getProcessor().setPixels(outPixels);

        return image;
    }

    /**
     *
     * @param rot
     * @param tilt
     */
    public ArrayList<String> getFiles(double rot, double tilt) {
        Point2d point = new Point2d(rot, tilt);
        Geometry.fixPoint(point);

        // File names are stored at startup so each cell contains a list of
        // the ones affecting it.
        SphereCell cell = getTableCell(point);

        return cell != null ? cell.fileNames : null;
    }

    private static void writeSelFile(String selFileName, ArrayList<String> fileNames) {
        MetaData md = new MetaData();

        md.addLabel(MDLabel.MDL_IMAGE);
        md.addLabel(MDLabel.MDL_ENABLED);

        for (int i = 0; i < fileNames.size(); i++) {
            long id = md.addObject();
            md.setValueString(MDLabel.MDL_IMAGE, fileNames.get(i), id);
            md.setValueInt(MDLabel.MDL_ENABLED, 1, id);
        }

        try {
            md.write(selFileName);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }
    }

    public ImageGeneric getProjection(ImageGeneric xmippVolume, double matrix[]) throws Exception {
        ImageGeneric projection = new ImageGeneric(ImageGeneric.Double);
        projection.resize(xmippVolume.getYDim(), xmippVolume.getXDim());

        IJ.showStatus(LABELS.MESSAGE_RETRIEVING_PROJECTION);
        Projection.projectVolume(xmippVolume, projection, matrix);
        //        projection.convert2Datatype(ImageGeneric_.SChar);

        return projection;
    }

    public String analyzeProjection(ImageGeneric xmippVolume, double matrix[], int w, int h) {
        String path = tempDir.getAbsolutePath() + File.separator;
        String projectionFileName = path + "projection.xmp";
        String selFileName = path + "analize.sel";
        String outputFile = path + "score.xmd";

        try {
            // * Gets projection.
            IJ.showStatus(LABELS.MESSAGE_RETRIEVING_PROJECTION);
            ImageGeneric projection = getProjection(xmippVolume, matrix);

            // * Writes projection to disk.
            projection.write(projectionFileName);

            // * Gets file names related to ROT and TILT
            // @TODO: getAngles...
            double rot = 0;
            double tilt = 0;
            ArrayList<String> files = getFiles(rot, tilt);

            // * Generate .sel file
            writeSelFile(selFileName, files);

            // * Calls xmipp_classify_...
            IJ.showStatus(LABELS.MESSAGE_ANALYZING_PROJECTION);

            String args = "-i " + selFileName
                    + " --ref " + projectionFileName
                    + //"--produceAligned", alignedFile,
                    " -o " + outputFile;

            System.out.println(" *** command: " + XMIPP_CLASSIFY_ANALYZE_CLUSTER + " " + args);
            //int result = Program.runByName(XMIPP_CLASSIFY_ANALYZE_CLUSTER, args);

            //System.out.println(" *** result: " + result);
            xmipp_classify_analyze_cluster(selFileName, projectionFileName, outputFile);
        } catch (Exception ex) {
            IJ.error("Error writing projection to temporary file: " + ex.getMessage());
            throw new RuntimeException(ex);
        }

        return outputFile;
    }

    // Example:
    // xmipp_classify_analyze_cluster
    //      -i xmipp_sel_file.sel
    //      --ref projection.xmp
    //      -o score.xmd
    private static void xmipp_classify_analyze_cluster(String selfile, String reffile, String outputFile) throws Exception {
        String command[] = new String[]{
            XMIPP_CLASSIFY_ANALYZE_CLUSTER,
            "-i", selfile,
            "--ref", reffile,
            //"--produceAligned", alignedFile,
            "-o", outputFile};

        System.out.print(" >>> ");
        for (int i = 0; i < command.length; i++) {
            System.out.print(command[i] + " ");
        }
        System.out.println();

        runCommand(command);
    }

    private static void runCommand(String command[]) throws Exception {
        //xmipp_classify_analyze_cluster -i xmipp_sel_file.sel -ref projection.xmp -produceAligned -oext ali.xmp -odir /tmp/projectionexplorer
        Process p = Runtime.getRuntime().exec(command);

        BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

        // read the output from the command
        String s;
        while ((s = stdInput.readLine()) != null) {
            System.out.println(s);
        }
    }
}
