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
	public long[] ids;
	public String[] mdBlocks;
	public String selectedBlock;
	// The following is only used in VolumeGallery mode
	public String selectedVol = "";
	public String[] volumes = null;

	public ArrayList<ColumnInfo> labels = null;
	//First label that can be rendered
	ColumnInfo ciFirstRender = null;
	public int zoom;
	public String filename;
//	public boolean galleryMode = true; // if false, is table model
//	public boolean volumeMode = false;
	public static final int MODE_GALLERY_MD = 1;
	public static final int MODE_GALLERY_VOL = 2;
	public static final int MODE_TABLE_MD = 3;
	
	//define min and max render dimensions
	public static final int MIN_SIZE = 16;
	public static final int MAX_SIZE = 256;
	
	//max dimension allowed to render images
	
	private int mode;
	public boolean showLabel = false;
	public boolean globalRender = false;
	public Param parameters;
	private int numberOfVols = 0;
	
	//flag to perform global normalization
	public boolean normalize = false;
	//flag to use geometry info
	public boolean useGeo = true;

	/**
	 * The constructor receive the filename of a metadata The metadata can also
	 * be passed, if null, it will be readed from filename
	 */
	public GalleryData(String fn, Param param, MetaData md) {
		try {
			selectedBlock = "";
			parameters = param;
			zoom = param.zoom;
			mode = MODE_GALLERY_MD;
			if (param.mode.equalsIgnoreCase(Param.OPENING_MODE_METADATA))
					mode = MODE_TABLE_MD;

			filename = fn;
			if (Filename.hasPrefix(fn)) {
				if (Filename.isMetadata(fn)) {
					selectedBlock = Filename.getPrefix(fn); // FIXME: validate block exists
					filename = Filename.getFilename(fn);
				}
			} 
			mdBlocks = MetaData.getBlocksInMetaDataFile(filename);

			if (mdBlocks.length > 1 && selectedBlock.isEmpty())
				selectedBlock = mdBlocks[0];

			if (md == null) {
				this.md = new MetaData();
				readMetadata(fn);
			} else {
				this.md = md;
				loadMd();
			}

		} catch (Exception e) {
			e.printStackTrace();
			md = null;
		}

	}// constructor GalleryData

	/** Return the name of the selected md block */
	public String getMdFilename() {
		if (selectedBlock.isEmpty())
			return filename;
		return String.format("%s@%s", selectedBlock, filename);
	}// function getMdFilename

	/** Load contents from a metadata already read */
	private void loadMd() throws Exception {
		ids = md.findObjects();
		loadLabels();
		numberOfVols = 0;
		volumes = null;
		if (isGalleryMode())
			mode = MODE_GALLERY_MD;
		
		if (hasRenderLabel()) {
			String imageFn = md.getValueString(ciFirstRender.getLabel(), md.firstObject());
			ImageGeneric image = new ImageGeneric(imageFn);
			
			//if (zoom == 0){ //default value
				int xdim = image.getXDim();
				int x = Math.min(Math.max(xdim, MIN_SIZE), MAX_SIZE);
				float scale = (float)x / xdim;
				zoom = (int) Math.ceil(scale * 100);
				DEBUG.printMessage(String.format("xdim: %d, x: %d, scale: %f, zoom: %d", xdim, x, scale, zoom));
			//}
				
			if (image.isVolume()) { // We are assuming all are volumes
									// or images, dont mix it
				if (isGalleryMode())
					mode = MODE_GALLERY_VOL;
				numberOfVols = md.size();
				volumes = new String[numberOfVols];
				for (int i = 0; i < numberOfVols; ++i)
					volumes[i] = md.getValueString(ciFirstRender.getLabel(), ids[i]);
				if (selectedVol.isEmpty())
					selectedVol = volumes[0];
			}
			image.destroy();
		}
		else //force this mode when there aren't render label
			mode = MODE_TABLE_MD;
		
	}// function loadMd
	
	/** Load labels info in md, 
	 * try to keep previous settings of render and visible
	 * on same columns */
	public void loadLabels(){
		ColumnInfo ci;
		try {
		int [] lab = md.getActiveLabels();
		ArrayList<ColumnInfo> newLabels = new ArrayList<ColumnInfo>(lab.length);
		ciFirstRender = null;
		ColumnInfo ciFirstRenderVisible = null;
		
		for (int i = 0; i < lab.length; ++i) {
			ci = new ColumnInfo(lab[i]);
			if (labels != null)
				for (ColumnInfo ci2: labels)
					if (ci.label == ci2.label){
						ci.updateInfo(ci2);
					}
			newLabels.add(ci);
			if (ciFirstRender == null && ci.allowRender)
				ciFirstRender = ci;
			if (ciFirstRenderVisible == null && ci.allowRender && ci.visible)
				ciFirstRenderVisible = ci;			
		}
		if (ciFirstRenderVisible != null)
			ciFirstRender = ciFirstRenderVisible;
		
		labels = newLabels;
		} catch (Exception e){
			e.printStackTrace();
		}
	}//function loadLabels

	/** Read metadata and store ids */
	private void readMetadata(String fn) {
		try {
			md.read(fn);
			loadMd();

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			md = null;
			ids = null;
		}
	}

	/** Select one of the blocks */
	public void selectBlock(String block) {
		selectedBlock = block; // FIXME: validate block exists
		readMetadata(getMdFilename());
	}

	public ImageGallery createModel() {
		try {
			switch (mode){
			case MODE_GALLERY_VOL:
				return new VolumeGallery(this);
			case MODE_GALLERY_MD:
				if (hasRenderLabel())
					return new MetadataGallery(this);
				//else fall in the next case
			case MODE_TABLE_MD:
				mode = MODE_TABLE_MD; //this is necessary when coming from previous case
				return new MetadataTable(this);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;

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

	public int getNumberOfVols() {
		return numberOfVols;
	}
	
	/** Return the mode of the gallery */
	public int getMode(){
		return mode;
	}
	
	/** Return true if there is a renderizable label in the metadata */
	public boolean hasRenderLabel(){
		return ciFirstRender != null;
	}
	
	//some mode shortcuts
	public boolean isGalleryMode(){return mode == MODE_GALLERY_MD || mode == MODE_GALLERY_VOL; }
	public boolean isVolumeMode() {return mode == MODE_GALLERY_VOL; }
	public boolean isTableMode() { return mode ==  MODE_TABLE_MD; }
	
	//utility function to change of mode
	public void changeMode() {
		if (isGalleryMode())
			mode = MODE_TABLE_MD;
		else if (numberOfVols > 0)
			mode = MODE_GALLERY_VOL;
		else
			mode = MODE_GALLERY_MD;
	}

	/** following function only should be used in VolumeGallery mode */
	public String getVolumeAt(int index) {
		return volumes[index];
	}

	public void selectVolume(String vol) {
		selectedVol = vol; // FIXME: Check it is valid
	}
	
	//Check if the underlying data has geometrical information
	public boolean containsGeometryInfo(){
		try {
			return md.containsGeometryInfo();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
	}

}// class GalleryData