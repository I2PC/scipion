import ij.ImagePlus;
import ij.process.FloatProcessor;
import xmipp.jni.ImageDouble;

/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

public class Converter {
	public static ImagePlus convertToImagePlus(ImageDouble img) {
		double [] imageData = img.getData();
		if (imageData == null)
			return null;
		
		// Create image
		int width=img.getXsize(), height=img.getYsize();
		FloatProcessor convertedImageProcessor = new FloatProcessor(width,height);
		int i=0;
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)	
					convertedImageProcessor.setf(x, y, (float) imageData[i++]);

		// normalize the image - done by default when calling FloatProcessor(array), here is our responsibility
		convertedImageProcessor.findMinAndMax();
		
		ImagePlus imagePlus = new ImagePlus("ImageDouble", convertedImageProcessor);
		return imagePlus;
	}
	
	public static ImageDouble convertToImageDouble(ImagePlus img) {
		if(img == null){
			Logger.debug("Null image");
			return null;
		}
		// Creates image 
		int width=img.getWidth(), height=img.getHeight(), numberOfProjections=img.getStackSize();
		ImageDouble imageDouble=new ImageDouble();

		int imageSize=width*height;
		double data[]=new double[imageSize*numberOfProjections];
		int i=0;
		// ImagePlus does not return the whole stack as a 3D array. 
		// Maybe this loop will be faster if block operations are used - like replacing the inner loop with a memcpy
		// Or, storing the float pixels array directly into the ImageDouble, with a new method like setProjection(float [])
		for(int p=1;p<=numberOfProjections;p++){
			FloatProcessor projection=(FloatProcessor) img.getStack().getProcessor(p);
			float [] pixels =(float[]) projection.getPixels();
			for (int j = 0; j < imageSize; j++)
					data[i++]=pixels[j]; 
		}
		
		try{
			// for tomograms we need N instead of Z, so Z = 1
			imageDouble.setData(width, height, 1, numberOfProjections, data);
		}catch (Exception ex){
			Logger.debug("convert ImagePlus->ImageDouble", ex);
		}
		return imageDouble;
	}
}
