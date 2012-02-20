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
	
	public ArrayList<ColumnInfo> labels;
	public int zoom;
	public String filename;
	public boolean galleryMode = true; // if false, is table model
	public boolean volumeMode = false; 
	public boolean showLabel = false;
	public Param parameters;
	private int numberOfVols = 0;
	

	/** The constructor receive the filename of a metadata 
	 * The metadata can also be passed, if null, it will be readed from filename
	 */
	public GalleryData(String fn, Param param, MetaData md) {
		try {
			selectedBlock = "";	
			parameters = param;
			zoom = param.zoom;
			galleryMode = param.mode.equalsIgnoreCase(Param.OPENING_MODE_GALLERY);
			
			//this should be moved to other place
			if (fn.contains("@")){
				String[] parts = fn.split("@");
				selectedBlock = parts[0]; //FIXME: validate block exists
				filename = parts[1];
			}
			else 
				filename = fn;
			mdBlocks = MetaData.getBlocksInMetaDataFile(filename);
			
			if (mdBlocks.length > 1 && selectedBlock.isEmpty())
				selectedBlock = mdBlocks[0];

			if (md == null){
				this.md = new MetaData();
				readMetadata(fn);
			}
			else {
				this.md = md;
				loadMd();
			}			
			
		} catch (Exception e) {
			e.printStackTrace();
			md = null;
		}

	}// constructor GalleryData
	
	/** Return the name of the selected md block */
	public String getMdFilename(){
		if (selectedBlock.isEmpty())
			return filename;
		return String.format("%s@%s", selectedBlock, filename);
	}//function getMdFilename
	
	/** Load contents from a metadata already read */
	private void loadMd() throws Exception{
		ids = md.findObjects();
	    volumeMode = false;
	    //if (galleryMode)
		if (md.containsLabel(MDLabel.MDL_IMAGE)) {
			String imageFn = md.getValueString(MDLabel.MDL_IMAGE,
					md.firstObject());
			ImageGeneric image = new ImageGeneric(imageFn);
			if (image.isVolume()) { // We are assuming all are volumes
									// or images, dont mix it
				volumeMode = true;
				numberOfVols = md.size();
				volumes = new String[numberOfVols];
				for (int i = 0; i < numberOfVols; ++i)
					volumes[i] = md.getValueString(MDLabel.MDL_IMAGE, ids[i]);
				if (selectedVol.isEmpty())
					selectedVol = volumes[0];
			}
			image.destroy();
		}
	    
	    if (!volumeMode){
	    	numberOfVols = 0;
	    	volumes = null;
	    }
		labels = ColumnInfo.createListFromMd(md);
	}//function loadMd
	
	/** Read metadata and store ids */
	private void readMetadata(String fn){
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
	
	/** Select one of the blocks*/
	public void selectBlock(String block){
		selectedBlock = block; //FIXME: validate block exists
		readMetadata(getMdFilename());
	}

	public ImageGallery createModel() {
		ImageGallery gallery = null;
		try {
			if (galleryMode) {
				gallery = (volumeMode ? new VolumeGallery(this)
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
	
	public int getNumberOfVols(){
		return numberOfVols;
	}
	
	/** following function only should be used in VolumeGallery mode */
	public String getVolumeAt(int index){
		return volumes[index];
	}
	
	public void selectVolume(String vol){
		selectedVol = vol; //FIXME: Check it is valid
	}
}// class GalleryData