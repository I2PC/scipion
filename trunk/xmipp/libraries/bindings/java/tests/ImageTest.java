package tests;

import xmipp.ImageGeneric;
import xmippij.XmippImageConverter;

public class ImageTest {

	public static void main(String args[]) {
		if (args.length < 1) {
			System.out.println("Usage: java ImageTest <xmipp_image_file>");
			System.exit(0);
		}

		try {
			String filename = args[0];

			// Reads image.
			System.out.println(" *** Reading image: " + filename);
			XmippImageConverter.convertToImagej(filename).show();
//			ImageGeneric image = new ImageGeneric(filename);
//
//			float pixels[] = image.getArrayFloat();
//			System.out.println(" *** Reading image (datatype: "
//					+ image.dataType + "):" + " x=" + image.xSize + " y="
//					+ image.ySize + " z=" + image.zSize + " n=" + image.nSize
//					+ " > datalength: " + pixels.length);
//
//			for (int j = 0; j < image.ySize && j < 3; j++) {
//				for (int i = 0; i < image.xSize && i < 3; i++) {
//					System.out.print(pixels[j * image.ySize + i] + " ");
//				}
//				System.out.println("...");
//			}
//			System.out.println("...");
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}
