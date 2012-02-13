package xmipp.viewer;

import ij.ImagePlus;
import xmipp.ij.XmippImageConverter;
import xmipp.jni.ImageGeneric;

public class VolumeGallery extends ImageGallery {
	public VolumeGallery(String fn, int zoom) throws Exception {
		super(fn, zoom);
		normalize = true; //volumes are displayed with global normalization by default
		calculateMinAndMax();
	}

	private static final long serialVersionUID = 1L;

	// Load initial dimensions
	protected ImageDimension loadDimension() throws Exception{
		ImageGeneric image = new ImageGeneric(filename); // read image header
		ImageDimension dim = new ImageDimension(image);
		image.destroy();
		return dim;
	}
	
	@Override
	protected double[] getMinAndMax() {
        try {
            ImageGeneric image = new ImageGeneric(filename);
            image.read(ImageGeneric.FIRST_IMAGE);
            double[] stats = image.getStatistics();
            image.destroy();
            return stats;
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return null;
	}

	@Override
	protected String getItemKey(int index) {
		try {
			return String.format("%d_(%d,%d)", index, thumb_width,
					thumb_height);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	@Override
	protected ImageItem createItem(int index, String key) {
		try {
			ImageGeneric image = new ImageGeneric(filename);
			ImagePlus imp = XmippImageConverter.convertToImageJ(image, thumb_width, thumb_height, index + 1, ImageGeneric.FIRST_IMAGE);
			image.destroy();
			return new ImageItem(key, key, imp);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

}
