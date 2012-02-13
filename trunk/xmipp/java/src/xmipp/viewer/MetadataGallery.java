package xmipp.viewer;

import ij.ImagePlus;
import xmipp.ij.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;

public class MetadataGallery extends ImageGallery {

	private static final long serialVersionUID = 1L;

	// MetaData with elements to display
	protected MetaData md;
	// Label to be rendered
	protected int renderLabel;
	protected int displayLabel;

	// Ids of objects stored in metadata
	long[] ids;

	public MetadataGallery(String fn, int zoom) throws Exception {
		super(fn, zoom);
	}

	// Load initial dimensions
	protected ImageDimension loadDimension() throws Exception {
		md = new MetaData(filename);
		ids = md.findObjects();
		renderLabel = MDLabel.MDL_IMAGE;
		displayLabel = MDLabel.MDL_IMAGE;
		ImageGeneric image = getImage(0);
		ImageDimension dim = new ImageDimension(image);
		// TODO: check this well, now asuming not volumes in metadata
		dim.setZDim(ids.length);
		image.destroy();
		return dim;
	}

	@Override
	protected ImageItem createItem(int index, String key) throws Exception {
		ImageGeneric image = getImage(index);
		image.readApplyGeo(md, ids[index], thumb_width, thumb_height);
		ImagePlus imp = XmippImageConverter.convertImageGenericToImageJ(image);
		String label = md.getValueString(renderLabel, ids[index]);
		return new ImageItem(key, label, imp);
	}

	@Override
	protected String getItemKey(int index) throws Exception {
		String filename = md.getValueString(renderLabel, ids[index]);
		return String.format("%s(%d,%d)", filename, thumb_width, thumb_height);
	}
	
	@Override
	public String getTitle() {
		return String.format("Metadata: %s (%d)", filename, n);
	}

	/**
	 * Create and returns the image generic from a given an index
	 * 
	 * @param index
	 *            Index of the image
	 * @return ImageGeneric created
	 * @throws Exception
	 *             if can not load image
	 */
	protected ImageGeneric getImage(int index) throws Exception {
		String imgFn = md.getValueString(renderLabel, ids[index]);
		return new ImageGeneric(imgFn);
	}

	@Override
	protected double[] getMinAndMax() {
		try {
			return md.getStatistics(false);
		} catch (Exception ex) {
			DEBUG.printException(ex);
		}
		return null;
	}
}
