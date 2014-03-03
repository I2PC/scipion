package xmipp.ij.commons;

import java.awt.Image;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import ij.IJ;
import ij.ImagePlus;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;


public class XmippUtil {

	private static XmippImageJ xij;

	public static XmippImageJ showImageJ(Tool tool) {
		if (IJ.getInstance() == null) {
			xij = new XmippImageJ();
			IJ.run("Install...",
					"install="
							+ Filename
									.getXmippPath("java/src/xmipp/ij/commons/XmippMacros.txt"));
		} else if (!xij.isVisible())
			xij.setVisible(true);
		return xij;
	}

	public static XmippImageJ getXmippImageJ() {
		return xij;
	}
        
        

	public static ImagePlus getImagePlus(String file) {
		try {

			ImageGeneric ig = new ImageGeneric(file);
			ig.read(ImageGeneric.FIRST_IMAGE);
			ImagePlus imp = XmippImageConverter.convertToImagePlus(ig);
			ig.destroy();

			return imp;
		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	public static Icon getImageIcon(ImagePlus imp, int width, int height)
	{

		Image image = imp.getImage().getScaledInstance(width, height, Image.SCALE_SMOOTH);
		Icon icon = new ImageIcon(image);

		return icon;
	}
        
        public static void executeCommand(String command) throws Exception {

            System.out.println(command);
            StringBuffer output = new StringBuffer();

            Process p;

            p = Runtime.getRuntime().exec(command);
            p.waitFor();
            BufferedReader reader
                    = new BufferedReader(new InputStreamReader(p.getInputStream()));

            output.append("Output\n");
            String line = "";
            while ((line = reader.readLine()) != null) {
                output.append(line + "\n");
            }
            reader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            output.append("Error\n");
            while ((line = reader.readLine()) != null) {
                output.append(line + "\n");
            }
            System.out.println(output.toString());
        }


}
