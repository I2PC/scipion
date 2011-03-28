package tests;

import xmipp.ImageDouble;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class PreviewTest {

    public static void main(String args[]) {
        if (args.length < 1) {
            System.out.println("Missing argument: A file is needed.");
	    System.exit(0);
        }

        try {
            String filename = args[0];
            ImageDouble image = new ImageDouble();
            image.readHeader(filename);
            image.printShape();

            int w = 60;//image.getXsize();
            int h = 60;//image.getYsize();
            int d = image.getZsize();
            long n = image.getNsize();

            //ImageStack is = null;
            image = new ImageDouble();

            for (int i = ImageDouble.FIRST_IMAGE; i <= n; i++) {
                for (int j = ImageDouble.FIRST_SLICE; j <= d; j++) {
                    image.readPreview(filename, w, h, j, i);

/*                    ImagePlus slice = ImageConverter.convertToImagej(image, filename);

                    if (is == null) {
                        is = new ImageStack(slice.getWidth(), slice.getHeight());
                    }

                    is.addSlice("", slice.getProcessor());*/
                    System.out.println(" @ " + i + " / " + j);// + ": " + slice.getProcessor().getPixelCount());
                }
            }
/*
            ImagePlus ip = new ImagePlus("", is);
            ip.show();*/
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }
}

