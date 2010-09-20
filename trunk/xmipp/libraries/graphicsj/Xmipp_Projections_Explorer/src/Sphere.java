
import constants.LABELS;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.Vector;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Sphere {

    private final static int R = 50;   // Sphere radius.
    private final static int ROT = Geometry.MAX_ROT;   // [0 - 359] = 360 elements.
    private final static int TILT = Geometry.MAX_TILT + 1;   // [0 - 180] = 181 elements.
    private final static int CRUST_THICKNESS = 2;
    private final static double MIN_SPREAD_VALUE = 1;   // Minimum value to spread after a hit.
    public final static int VALUE = 255;   // Initial value for hits.
    private SphereCell sphereMap[][];
    private final static String TEMPDIR_PATH = System.getProperty("java.io.tmpdir") + "/explorer";
    private File tempDir;
    //private final static String work_dir = "unknown! - @TODO: directorio actual.";
    private final static String XMIPP_CLASSIFY_ANALYZE_CLUSTER = "xmipp_classify_analyze_cluster";

    public Sphere() {
        /*        sphereMap = new SphereCell[ROT][TILT];
        for (int i = 0; i < sphereMap.length; i++) {
        for (int j = 0; j < sphereMap[0].length; j++) {
        sphereMap[i][j] = new SphereCell();
        }
        }*/
    }

    public Sphere(String fileName) {
        sphereMap = new SphereCell[ROT][TILT];
        for (int i = 0; i < sphereMap.length; i++) {
            for (int j = 0; j < sphereMap[0].length; j++) {
                sphereMap[i][j] = new SphereCell();
            }
        }

        load(fileName);

        equalizeValues(sphereMap);

        // Creates temporary work directory.
        tempDir = new File(TEMPDIR_PATH);
        try {
            Runtime.getRuntime().exec("rm -rf " + TEMPDIR_PATH);
            Runtime.getRuntime().exec("mkdir " + TEMPDIR_PATH);
        } catch (IOException ioex) {
            IJ.write(ioex.getMessage());
            ioex.printStackTrace();
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

    private void load(String filename) {
        try {
            IJ.showStatus(LABELS.MESSAGE_LOADING_EULER_ANGLES_FILE);
            File file = new File(filename);
            BufferedReader br = new BufferedReader(new FileReader(file));

            int fileSize = (int) file.length();
            int currentSize = 0;

            String baseDir = file.getParent() + File.separator;

            // Skips header.
            br.readLine();
            br.readLine();
            br.readLine();

            String line;
            String imageFileName;

            while ((line = br.readLine()) != null) {    // Reads line.
                StringTokenizer input = new StringTokenizer(line, " ");

                imageFileName = input.nextToken();  // Parses file name.

                input.nextToken();  // Parses enabled.

                // Parses ROT and TILT
                double rot = Double.parseDouble(input.nextToken());
                double tilt = Double.parseDouble(input.nextToken());

                // Builds absolute path (if it's not absolute already)
                if (!imageFileName.startsWith("/")) {
                    imageFileName = baseDir + imageFileName;
                }

                double angles[] = Geometry.normalizeAngles(rot, tilt);
                Point2d hit = new Point2d(angles[0], angles[1]);

                spreadHit(hit, imageFileName);

                // Calcualtes and shows progress...
                currentSize += imageFileName.getBytes().length;
                IJ.showProgress(currentSize, fileSize);
            }
        } catch (Exception ex) {
            IJ.write(ex.getMessage());
            ex.printStackTrace();
        }
    }

    private void spreadHit(Point2d hitted, String imageFileName) {
        setLocked(hitted, true);   // Locks hitted point...

        double value = addHit(hitted, hitted, imageFileName);    // ...adds a hit...

        Vector<Point2d> points = new Vector<Point2d>();
        spreadHit(hitted, hitted, value, imageFileName, points);    // ...and spreads it...

        setLocked(hitted, false);   // Finally unlocks all the points.
        for (int i = 0; i < points.size(); i++) {
            setLocked(points.elementAt(i), false);
        }
    }

    private void spreadHit(Point2d hitted, Point2d current, double value, String imageFileName, Vector<Point2d> points) {
        if (value > MIN_SPREAD_VALUE) {
            // Calculates all neighbours.
            Vector<Point2d> neighbours = getNeighbours(current);
            double values[] = new double[neighbours.size()];    // Individual values.

            // For each neighbour...
            for (int i = 0; i < neighbours.size(); i++) {
                Point2d neighbour = neighbours.elementAt(i);

                // ...adds value (stores it to spread it later)...
                values[i] = addHit(neighbour, hitted, imageFileName);

                // ...and adds point to the locked list.
                points.add(neighbour);
            }

            // After every neighbour is locked an calculated, spreads it
            // recursively (keeping the last value added to know when to stop it.)
            for (int i = 0; i < neighbours.size(); i++) {
                spreadHit(hitted, neighbours.elementAt(i), values[i], imageFileName, points);
            }
        }
    }

    /**
     * Gets all (not locked) neighbour points.
     * @param center
     * @return
     */
    private Vector<Point2d> getNeighbours(Point2d center) {
        Vector<Point2d> points = new Vector<Point2d>();

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

    public ImagePlus build() {
        int W, H, D;
        W = H = D = 2 * R + 1;  // For each dimension: r - 0 - r = 2r + 1

        ImageProcessor ip = new ColorProcessor(W, H);
        ImagePlus image = new ImagePlus("Sphere", ip);

        // Room for output ~ [slice][image]
        int[][] stack = new int[D][W * H];

        buildSphere(stack, W, H, D, R, CRUST_THICKNESS);

        // Builds the stack.
        ImageStack outputStack = new ImageStack(W, H, D);
        for (int slice = 0; slice < stack.length; slice++) {
            outputStack.setPixels(stack[slice], slice + 1);
        }

        image.setStack("Sphere", outputStack);

        return image;
    }

    private void buildSphere(int stack[][], int W, int H, int D, int r, int thickness) {
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
    public Vector<String> getFiles(double rot, double tilt) {
        Point2d point = new Point2d(rot, tilt);
        Geometry.fixPoint(point);

        // File names are stored at startup so each cell contains a list of
        // the ones affecting it.
        SphereCell cell = getTableCell(point);

        return cell != null ? cell.fileNames : null;
    }

    private static void writeSelFile(String selFileName, Vector<String> fileNames) {
        try {
            File selFile = new File(selFileName);
            //selFile.deleteOnExit();

            BufferedWriter out = new BufferedWriter(new FileWriter(selFile, false));

            // Writes header
            out.write("; XMIPP_3 * column_format * \n");
            out.write(";\n");
            out.write("; image enabled \n");

            // Writes files.
            for (int i = 0; i < fileNames.size(); i++) {
                out.write(fileNames.elementAt(i) + "\t1\n");
            }

            out.close();
        } catch (IOException e) {
            IJ.write("Cannot write temporary file to: " + System.getProperty("java.io.tmpdir"));
            IJ.write(e.getMessage());
        }
    }

    public ImagePlus getProjection(MultidimArrayd xmippVolume, double rot, double tilt, int w, int h) {
        ImageDouble id = new ImageDouble(w, h);

        Projection p = new Projection();
        p.reset(h, w);

        IJ.showStatus(LABELS.MESSAGE_RETRIEVING_PROJECTION);
        XmippData.project_Volume(xmippVolume, p, h, w, rot, tilt, 0);

        ImagePlus projection = ImageConverter.xmipp2imagej(p, w, h);   // To imageJ.

        p.delete();

        return projection;
    }

    public String[] analyzeProjection(MultidimArrayd xmippVolume, double rot, double tilt, int w, int h) {
        String pcaFile = "";
        String outliersFile = "";

        try {
            // * Gets projection.
            Projection p = new Projection();
            p.reset(h, w);

            IJ.showStatus(LABELS.MESSAGE_RETRIEVING_PROJECTION);
            XmippData.project_Volume(xmippVolume, p, h, w, rot, tilt, 0);

            // * Write projection to disk.
            String path = tempDir.getAbsolutePath() + File.separator;
            String projectionFileName = path + "projection";
            String projectionExt = ".xmp";
            String selFileName = path + "xmipp_sel_file.sel";

            FileName fn = new FileName(projectionFileName, projectionExt);
            p.write(fn);

            p.delete();

            // * Gets file names related to ROT and TILT
            Vector<String> files = getFiles(rot, tilt);

            // * Generate .sel file
            //IJ.write("File: " + selFileName);
            writeSelFile(selFileName, files);

            // * Call xmipp_classify_...
            IJ.showStatus(LABELS.MESSAGE_ANALYZING_PROJECTION);
            xmipp_classify_analyze_cluster(selFileName, projectionFileName + projectionExt, "ali.xmp");

            // * Builds "good" and "bad" images so aligned images can be shown at table.
            pcaFile = path + "xmipp_sel_file_pca.sel";
            outliersFile = path + "xmipp_sel_file_outliers.sel";
        } catch (Exception ioex) {
            IJ.write("Error writing projection to temporary file. ");
            IJ.write(ioex.getMessage());
            ioex.printStackTrace();
        }

        return new String[]{pcaFile, outliersFile};
    }

    private static void xmipp_classify_analyze_cluster(String selfile, String reffile, String oext) {
        try {
            //xmipp_classify_analyze_cluster -i xmipp_sel_file.sel -ref projection.xmp -produceAligned -oext ali.xmp -odir /tmp/projectionexplorer
            Process p = Runtime.getRuntime().exec(
                    "/home/juanjo/xmipp/bin/"
                    + XMIPP_CLASSIFY_ANALYZE_CLUSTER
                    + " -i " + selfile
                    + " -ref " + reffile
                    + " -produceAligned -oext " + oext
                    + " -odir " + TEMPDIR_PATH);

            BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

            // read the output from the command
            String s;
            while ((s = stdInput.readLine()) != null) {
                System.out.println(s);
            }
        } catch (IOException ioex) {
            IJ.write("Error executing " + XMIPP_CLASSIFY_ANALYZE_CLUSTER);
            IJ.write(ioex.getMessage());
        }
    }
}
