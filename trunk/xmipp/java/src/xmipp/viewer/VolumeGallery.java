package xmipp.viewer;

import ij.ImagePlus;
import xmipp.ij.XmippImageConverter;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;

public class VolumeGallery extends ImageGallery {
	protected String volFn;
	protected long volNumber;
	ImageGeneric volHeader;
	
	public VolumeGallery(GalleryData data) throws Exception {
		super(data);
		data.normalize = true; // volumes are displayed with global normalization by default
		volFn = Filename.getFilename(data.selectedVol);
		volNumber = Filename.getNimage(data.selectedVol);
		calculateMinAndMax();
	}

	private static final long serialVersionUID = 1L;

	// Load initial dimensions
	protected ImageDimension loadDimension() throws Exception {
		volHeader = new ImageGeneric(data.selectedVol); // read image header
		ImageDimension dim = new ImageDimension(volHeader);
		return dim;
	}

	@Override
	protected double[] getMinAndMax() {
		try {
			ImageGeneric image = new ImageGeneric(volFn);
			image.read(volNumber);
			double[] stats = image.getStatistics();
			image.destroy();
			return stats;
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return null;
	}

	@Override
	protected String getItemKey(int index) throws Exception {
		return String.format("%d_(%d,%d)", index, thumb_width, thumb_height);
	}
	
	@Override
	public String getTitle() {
		return String.format("Volume: %s (%d x %d x %d)", data.selectedVol, image_width, image_height, n);
	}

	@Override
	protected ImageItem createItem(int index, String key) throws Exception {
		ImagePlus imp = XmippImageConverter.readImageGenericToImageJ(volHeader, thumb_width,
				thumb_height, index + 1, volNumber);
		String label = String.format("%d", index + 1);
		return new ImageItem(key, label, imp);
	}

	@Override
	public ImagePlus getImagePlus() {
		try {
			return XmippImageConverter.readImageGenericToImageJ(volHeader);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

}
