package tests;

import xmipp.ImageDouble;

public class ImageTest {

    public static void main(String args[]) {
        if (args.length < 1) {
            System.out.println("Usage: java ImageTest <xmipp_image_file>");
            System.exit(0);
        }

        try {
            String file = args[0];

            // Reads image.
            ImageDouble image = new ImageDouble();
            image.readStack(file);
            System.out.println(" *** Full image: " + file);
            image.printShape();

            // Reads preview.
            ImageDouble preview = new ImageDouble();
            preview.readPreview(file, 80, 80);
            System.out.println(" *** Image preview: " + file);
            preview.printShape();

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
