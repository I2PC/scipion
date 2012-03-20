package xmipp.viewer.models;

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
	// First label that can be rendered
	ColumnInfo ciFirstRender = null;
	public int zoom;
	public String filename;
	// public boolean galleryMode = true; // if false, is table model
	// public boolean volumeMode = false;
	public enum Mode { GALLERY_MD, GALLERY_VOL, 
		                 TABLE_MD, GALLERY_ROTSPECTRA};
//	public static final int Mode.GALLERYGALLERY_MD = 1;
//	public static final int MODE_GALLERY_VOL = 2;
//	public static final int Mode.GALLERY_MD = 3;

	// define min and max render dimensions
	public static int MIN_SIZE = 16;
	public static int MAX_SIZE = 256;

	// max dimension allowed to render images

	private Mode mode;
	public boolean showLabel = false;
	public boolean globalRender;
	public Param parameters;
	private int numberOfVols = 0;

	// flag to perform global normalization
	public boolean normalize = false;
	// flag to use geometry info
	public boolean useGeo = true;
	// flag to wrapping
	public boolean wrap = true;
	// flag to check if is 2d classification
	public boolean is2dClassification = false;
	// Store the selection state for each item
	protected boolean[] selection;
	public int[] classes;

	/**
	 * The constructor receive the filename of a metadata The metadata can also
	 * be passed, if null, it will be readed from filename
	 */
	public GalleryData(String fn, Param param, MetaData md) {
		try {
			selectedBlock = "";
			parameters = param;
			zoom = param.zoom;
			globalRender = param.renderImages;
			mode = Mode.GALLERY_MD;
			
			if (param.mode.equalsIgnoreCase(Param.OPENING_MODE_METADATA))
				mode = Mode.TABLE_MD;
			else if (param.mode.equalsIgnoreCase(Param.OPENING_MODE_ROTSPECTRA))
				mode = Mode.GALLERY_ROTSPECTRA;

			filename = fn;
			if (Filename.hasPrefix(fn)) {
				if (Filename.isMetadata(fn)) {
					selectedBlock = Filename.getPrefix(fn); // FIXME: validate
															// block exists
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
		useGeo = containsGeometryInfo();
		selection = new boolean[ids.length];
		is2dClassification = checkifIs2DClassificationMd();
		
		if (is2DClassificationMd()){
			classes = new int[ids.length];
			for (int i = 0; i < ids.length; ++i)
				classes[i] = -1;
		}

		if (isRotSpectraMd())
		{
			mode = Mode.GALLERY_ROTSPECTRA;
			if (zoom == 0)
				zoom = 100;
			return;
		}
		
		if (isGalleryMode())
			mode = Mode.GALLERY_MD;

		if (hasRenderLabel()) {
			int renderLabel = ciFirstRender.getLabel();
			ImageGeneric image = null;
			String imageFn;
			for (int i = 0; i < ids.length && image == null; ++i) {
				imageFn = md.getValueString(renderLabel, ids[i]);
				if (Filename.exists(imageFn))
					image = new ImageGeneric(imageFn);
			}
			if (image != null) { // Image file was found to render
				if (zoom == 0) { // if default value, calculate zoom
					//If in micrograph mode, reduce the MAX_SIZE constant
					if (md.containsMicrographsInfo())
						MAX_SIZE /= 2;
					int xdim = image.getXDim();
					int x = Math.min(Math.max(xdim, MIN_SIZE), MAX_SIZE);
					float scale = (float) x / xdim;
					zoom = (int) Math.ceil(scale * 100);
				}

				if (image.isVolume()) { // We are assuming all are volumes
										// or images, dont mix it
					if (isGalleryMode())
						mode = Mode.GALLERY_VOL;
					numberOfVols = md.size();
					volumes = new String[numberOfVols];

					for (int i = 0; i < numberOfVols; ++i)
						volumes[i] = md.getValueString(
								ciFirstRender.getLabel(), ids[i]);
					if (selectedVol.isEmpty())
						selectVolume(volumes[0]);
				}
				image.destroy();
			}
		} else {
			// force this mode when there aren't render label
				mode = Mode.TABLE_MD;
		}

	}// function loadMd

	/**
	 * Load labels info in md, try to keep previous settings of render and
	 * visible on same columns
	 */
	public void loadLabels() {
		ColumnInfo ci;
		try {
			int[] lab = md.getActiveLabels();
			ArrayList<ColumnInfo> newLabels = new ArrayList<ColumnInfo>(
					lab.length);
			ciFirstRender = null;
			ColumnInfo ciFirstRenderVisible = null;

			for (int i = 0; i < lab.length; ++i) {
				ci = new ColumnInfo(lab[i]);
				if (labels != null) {
					for (ColumnInfo ci2 : labels)
						if (ci.label == ci2.label)
							ci.updateInfo(ci2);
				} else if (ci.allowRender)
					ci.render = globalRender;
				newLabels.add(ci);
				if (ciFirstRender == null && ci.allowRender)
					ciFirstRender = ci;
				if (ciFirstRenderVisible == null && ci.allowRender
						&& ci.visible)
					ciFirstRenderVisible = ci;
			}
			if (ciFirstRenderVisible != null) {
				ciFirstRender = ciFirstRenderVisible;
				// Add MDL_ENABLED if not present
				if (!md.containsLabel(MDLabel.MDL_ENABLED)) {
					newLabels.add(0, new ColumnInfo(MDLabel.MDL_ENABLED));
					md.addLabel(MDLabel.MDL_ENABLED);
					for (long id : ids)
						md.setEnabled(true, id);
				}
			}

			labels = newLabels;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}// function loadLabels

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
			switch (mode) {
			case GALLERY_VOL:
				return new VolumeGallery(this);
			case GALLERY_MD:
				if (hasRenderLabel())
					return new MetadataGallery(this);
				// else fall in the next case
			case TABLE_MD:
				mode = Mode.TABLE_MD; // this is necessary when coming from
				// previous case
				if (!md.isColumnFormat())
					return new MetadataRow(this);
				if (md.containsMicrographsInfo()) 
					return new MicrographsTable(this);
				return new MetadataTable(this);
			case GALLERY_ROTSPECTRA:
				return new RotSpectraGallery(this);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public int getNumberOfBlocks() {
		return mdBlocks.length;
	}

	public int getNumberOfVols() {
		return numberOfVols;
	}

	/** Return the mode of the gallery */
	public Mode getMode() {
		return mode;
	}

	/** Return true if there is a renderizable label in the metadata */
	public boolean hasRenderLabel() {
		return ciFirstRender != null;
	}

	/** Return the label that is used for rendering */
	public int getRenderLabel() {
		return ciFirstRender.getLabel();
	}
	
	/** Return true if the gallery mode is allowed */
	public boolean allowGallery(){
		return hasRenderLabel() || isRotSpectraMd();
	}

	// some mode shortcuts
	public boolean isGalleryMode() {
		return mode == Mode.GALLERY_MD || 
				mode == Mode.GALLERY_VOL ||
				mode == Mode.GALLERY_ROTSPECTRA;
	}

	public boolean isVolumeMode() {
		return mode == Mode.GALLERY_VOL;
	}

	public boolean isTableMode() {
		return mode == Mode.TABLE_MD;
	}
	
	public boolean isRotSpectraMode(){
		return mode == Mode.GALLERY_ROTSPECTRA;
	}
	
	public boolean isMicrographsMode(){
		return md.containsMicrographsInfo();
	}

	// utility function to change of mode
	public void changeMode() {
		if (isGalleryMode())
			mode = Mode.TABLE_MD;
		else if (isRotSpectraMd())
			mode = Mode.GALLERY_ROTSPECTRA;
		else if (numberOfVols > 0)
			mode = Mode.GALLERY_VOL;
		else
			mode = Mode.GALLERY_MD;
	}

	/** following function only should be used in VolumeGallery mode */
	public String getVolumeAt(int index) {
		return volumes[index];
	}

	public void selectVolume(String vol) {
		selectedVol = vol; // FIXME: Check it is valid
	}

	// Check if the underlying data has geometrical information
	public boolean containsGeometryInfo() {
		try {
			return md.containsGeometryInfo();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
	}

	/** Check if an item is enabled or not */
	public boolean isEnabled(int index) {
		try {
			if (!isVolumeMode()) // slices in a volume are always enabled
				return md.getEnabled(ids[index]);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return true;
	}

	/** Set enabled state */
	public void setEnabled(int index, boolean value) {
		try {
			if (!isVolumeMode()) // slices in a volume are always enabled
				md.setEnabled(value, ids[index]);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/** This is only needed for metadata table galleries */
	public boolean isFile(int col) {
		try {
			return MetaData.isPathField(labels.get(col).getLabel());
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}
	
	public boolean isImageFile(int col) {
		return isImageFile(labels.get(col));		
	}
	
	public boolean isImageFile(ColumnInfo ci) {
		try {
			return MetaData.isImage(ci.getLabel());
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}
	
	public boolean is2DClassificationMd(){
		return is2dClassification;
	}
	/** Return true if current metadata comes from 2d classification */
	public boolean checkifIs2DClassificationMd(){
		try {
			if (!selectedBlock.equalsIgnoreCase("classes"))
				return false;
			int n = md.size();
			for (long id: ids) {
				int ref = md.getValueInt(MDLabel.MDL_REF, id);
				String s = String.format("class%06d_images", ref);
				if (!containsBlock(s))
					return false;
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return true;
	}
	
	/** Get the assigned class of some element */
	public int getItemClass(int index){
		if (is2dClassification)
			return classes[index];
		return -1;
	}
	
	/** Set the class of an element */
	public void setItemClass(int index, int classNumber){
		classes[index] = classNumber;
	}
	
	public void setSelectionClass(int classNumber){
		
	}
	
	/** Return true if current metadata is a rotspectra classes */
	public boolean isRotSpectraMd(){
		String fnVectors = filename.replace("classes", "vectors");
		String fnVectorsData = fnVectors.replace(".xmd", ".vec");
		//TODO: CHECK if is a classification md
		if (is2DClassificationMd() &&
			Filename.exists(fnVectors) && 
			Filename.exists(fnVectorsData))
			return true;
		return false;
	}
	
	/** Check if a block is present, ignore case */
	public boolean containsBlock(String block){
		for (String b: mdBlocks)
			if (b.equalsIgnoreCase(block))
				return true;
		return false;
	}
	
	public String getValueFromCol(int index, int col){
		return getValueFromCol(index, labels.get(col));
	}
	
	public String getValueFromCol(int index, ColumnInfo ci){
		try {
			return md.getValueString(ci.getLabel(), ids[index]);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}	
	
	public String getValueFromLabel(int index, int label){
		try {
			return md.getValueString(label, ids[index]);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
}// class GalleryData