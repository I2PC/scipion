package xmipp.viewer;

import java.util.ArrayList;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;

/** This class will serve to store important data about the gallery */
public class GalleryData {
	public MetaData md;
	public String[] mdBlocks;
	public ArrayList<ColumnInfo> labels;
	public int zoom;
	public String filename;
	public boolean galleryMode = true; // if false, is table model
	public boolean showLabel = false;
	public int numberOfVols = 0;

	/** The constructor receive the filename of a metadata */
	public GalleryData(String filename, Param param) {
		try {
			mdBlocks = MetaData.getBlocksInMetaDataFile(filename);
			md = new MetaData(filename);
			labels = ColumnInfo.createListFromMd(md);
			zoom = param.zoom;
			this.filename = filename;
			galleryMode = param.mode
					.equalsIgnoreCase(Param.OPENING_MODE_GALLERY);
		} catch (Exception e) {
			e.printStackTrace();
			md = null;
		}

	}// constructor GalleryData

	public ImageGallery createModel() {
		ImageGallery gallery = null;
		try {
			if (galleryMode) {
				boolean volGallery = false;
				if (md.containsLabel(MDLabel.MDL_IMAGE)) {
					String fn = md.getValueString(MDLabel.MDL_IMAGE,
							md.firstObject());
					ImageGeneric image = new ImageGeneric(fn);
					if (image.isVolume()) { // We are assuming all are volumes
											// or images, dont mix it
						volGallery = true;
						numberOfVols = md.size();
					}
				}
				gallery = (volGallery ? new VolumeGallery(this)
						: new MetadataGallery(this));
			} else
				gallery = new MetadataTable(this);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return gallery;

		// if (Filename.isVolume(filename))
		// gallery = new VolumeGallery(data);
		// else {
		// DEBUG.printMessage(parameters.mode);
		// if (parameters.mode.equalsIgnoreCase(Param.OPENING_MODE_GALLERY))
		// {//if (Filename.isStack(filename))
		// DEBUG.printMessage("creatingMetadataGallery ");
		// gallery = new MetadataGallery(data);
		// }
		// //else (Filename.isMetadata(filename))
		// else if
		// (parameters.mode.equalsIgnoreCase(Param.OPENING_MODE_METADATA))
		//
		// }
	}

	public int getNumberOfBlocks() {
		return mdBlocks.length;
	}
}// class GalleryData