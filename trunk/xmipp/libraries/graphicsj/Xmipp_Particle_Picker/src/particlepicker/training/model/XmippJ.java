package particlepicker.training.model;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.util.Tools;

import java.io.File;
import java.util.List;

import xmipp.jni.ImageGeneric;
import xmippij.XmippImageConverter;

public class XmippJ {
	
	
	public static String saveTempImageStack(List<ImagePlus> imgs)
	{
		if(imgs.size() == 0)
			throw new IllegalArgumentException("No images provided");
			try {
				ImageStack stack = null;
				for (ImagePlus iplus :imgs) {
					if (stack == null)
						stack = new ImageStack(iplus.getWidth(), iplus.getHeight());
					stack.addSlice("", iplus.getProcessor().convertToFloat());
				}
				
				ImagePlus ipstack = new ImagePlus("", stack);
				
				//ImageGeneric idouble = new ImageGeneric();

//		        int w = ipstack.getWidth();
//		        int h = ipstack.getHeight();
//		        int d = ipstack.getStackSize();
	//
//		        double data[] = new double[w * h * d];
//		        for (int i = 0; i < d; i++) {
//		            float slice[] = (float[]) ipstack.getStack().getProcessor(i + 1).getPixels();
//		            System.arraycopy(Tools.toDouble(slice), 0, data, i * w * h, w * h);
//		        }
//		        try {
//		            image.setData(w, h, d, data);
//		        } catch (Exception ex) {
//		            ex.printStackTrace();
//		        }
//				
//				System.out.println(ipstack.getWidth());
//				System.out.println(ipstack.getHeight());
//				System.out.println(ipstack.getStackSize());
//				System.out.println(data.length);
//				System.out.println("-------------------");


				
				
				/*idouble.setData(ipstack.getWidth(), 
						ipstack.getHeight(), 
						ipstack.getStackSize(),
						Tools.toDouble((float[])ipstack.getProcessor().getPixels()));*/
				ImageGeneric idouble = XmippImageConverter.convertToXmipp(ipstack);
				File dir = new File(System.getProperty("user.dir"));
				
				File file = File.createTempFile("xmipp", ".stk", dir);
				idouble.write(file.getAbsolutePath());
				return file.getAbsolutePath();
			} catch (Exception e) {
				e.printStackTrace();
				IJ.error(e.getMessage());
			}
			return null;
		}


}
