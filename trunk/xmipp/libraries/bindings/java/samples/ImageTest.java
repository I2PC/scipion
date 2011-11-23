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
			//readArray(filename);
			setArray(filename);
			
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	public static void showImagej(String filename) throws Exception
	{
		// Reads image.
		System.out.println(" *** Reading image: " + filename);
		XmippImageConverter.convertToImagej(filename).show();
	}
	
	public static void readArray(String filename) throws Exception
	{
		ImageGeneric image = new ImageGeneric(filename);

		float pixels[] = image.getArrayFloat();
		System.out.println(" *** Reading image (datatype: "
				+ image.dataType + "):" + " x=" + image.xSize + " y="
				+ image.ySize + " z=" + image.zSize + " n=" + image.nSize
				+ " > datalength: " + pixels.length);
		System.out.println(" *** Array readed: ");
		for (int j = 0; j < image.ySize && j < 3; j++) {
			for (int i = 0; i < image.xSize && i < 3; i++) {
				System.out.print(pixels[j * image.ySize + i] + " ");
				
			}
			System.out.println("...");
		}
		System.out.println("...");
	}
	
	public static void setArray(String filename)throws Exception
	{
		System.out.println("Before constructor");
		ImageGeneric img = new ImageGeneric(filename);
		System.out.println("Reading filename: " + filename);
		img.readData(filename);
		System.out.println("Getting array...");
		float[] data = img.getArrayFloat();
		System.out.println("After getArray");
		for (int i = 0; i < img.xSize; i++)
			for (int j = 0; j < img.ySize; j++)
			{
				if (j < 16)
					data[i*img.xSize+j] = 0;
			}
		System.out.println("Before setArray");
		
		img.setArrayFloat(data);
		System.out.println("Before writing image");
		img.write("kk.xmp");
	}
}
