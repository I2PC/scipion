package tests;

import xmipp.ImageGeneric;
import xmipp.Projection;

public class ProjectionTest {

    public static void main(String args[]) {
        if (args.length < 2) {
            System.out.println("Usage: java ImageTest <xmipp_image_file> <output_projection_file>");
            System.exit(0);
        }

        try {
            String file = args[0];
            String outfile = args[1];

            // Loads volume.
            ImageGeneric image = new ImageGeneric();
            image.readData(file);
            image.setXmippOrigin();

            image.printShape();

            // Creates a projection and resets it.
            Projection projection = new Projection();
            projection.reset(image.xSize, image.ySize);

            // Retrieve projection for given angles.
            double rot = 45;
            double tilt = 90;
            double psi = 0;
            Projection.projectVolume(image, projection, rot, tilt, psi);

            // Save projection to a file.
            projection.write(outfile);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
