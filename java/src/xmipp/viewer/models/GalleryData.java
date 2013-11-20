/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

package xmipp.viewer.models;

import java.awt.Color;
import java.awt.Window;
import java.io.File;
import java.util.ArrayList;
import java.util.Date;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.jni.MDRow;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;
import xmipp.utils.XmippStringUtils;

/** This class will serve to store important data about the gallery */
public class GalleryData {
	public MetaData md;
	public long[] ids;
	public String[] mdBlocks = null;
	public String selectedBlock;
	// The following is only used in VolumeGallery mode
	public String selectedVolFn = "";
	public String commonVolPrefix = "";
	public String[] volumes = null;

	public ArrayList<ColumnInfo> labels = null;
	// First label that can be rendered
	ColumnInfo ciFirstRender = null;
	public int zoom;
	private String filename;
	public int resliceView;

	public enum Mode {
		GALLERY_MD, GALLERY_VOL, TABLE_MD, GALLERY_ROTSPECTRA
	};

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
	public boolean useGeo;
	// flag to wrapping
	public boolean wrap;
	// flag to check if is 2d classification
	public boolean is2dClassification = false;
	// Store the selection state for each item
	public boolean[] selection;
	// Array with all ClassInfo
	public ArrayList<ClassInfo> classesArray;
	// ClassInfo reference for each element
	public ClassInfo[] classes;
	// Flags to check if md or classes has changed
	private boolean hasMdChanges, hasClassesChanges;
	public Window window;

	/**
	 * The constructor receive the filename of a metadata The metadata can also
	 * be passed, if null, it will be readed from filename
	 * 
	 * @param jFrameGallery
	 */
	public GalleryData(Window window, String fn, Param param, MetaData md) {
		this.window = window;
		try {
			selectedBlock = "";
			parameters = param;
			zoom = param.zoom;
			globalRender = param.renderImages;
			mode = Mode.GALLERY_MD;
			resliceView = param.resliceView;
			useGeo = param.useGeo;
			wrap = param.wrap;

			if (param.mode.equalsIgnoreCase(Param.OPENING_MODE_METADATA))
				mode = Mode.TABLE_MD;
			else if (param.mode.equalsIgnoreCase(Param.OPENING_MODE_ROTSPECTRA))
				mode = Mode.GALLERY_ROTSPECTRA;

			setFileName(fn);

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

	public ArrayList<ColumnInfo> getColumns() {
		return labels;
	}

	public void setRenderColumn(ColumnInfo ci) {
		ciFirstRender = ci;
	}

	public ColumnInfo getRenderColumn() {
		return ciFirstRender;
	}

	/** Return the name of the selected md block */
	public String getMdFilename() {
		if (selectedBlock.isEmpty())
			return filename;
		return String.format("%s@%s", selectedBlock, filename);
	}// function getMdFilename

	public String getMdSaveFileName() {
		if (filename == null)
			return null;

		String savefn;
		if (selectedBlock.isEmpty())
			savefn = filename;
		else
			savefn = String.format("%s@%s", selectedBlock, filename);
		String ext;
		if (savefn.contains(".")) {
			ext = savefn.substring(savefn.lastIndexOf("."));
			if (ext.equals(".stk"))
				return savefn.replace(ext, ".xmd");
		} else
			savefn = savefn + ".xmd";
		return savefn;
	}

	public void setFileName(String file) {
		filename = file;

		if (file != null) {
			if (Filename.hasPrefix(file)) {
				if (Filename.isMetadata(file)) {
					selectedBlock = Filename.getPrefix(file); // FIXME:
																// validate
																// block exists
					filename = Filename.getFilename(file);
				}
			}
			mdBlocks = MetaData.getBlocksInMetaDataFile(filename);

			if (mdBlocks.length >= 1 && selectedBlock.isEmpty())
				selectedBlock = mdBlocks[0];
		}

	}

	/** Load contents from a metadata already read */
	public void loadMd() throws Exception {
		ids = md.findObjects();
		loadLabels();
		numberOfVols = 0;
		volumes = null;
		if (!containsGeometryInfo())
			useGeo = false;
		selection = new boolean[ids.length];
		is2dClassification = checkifIs2DClassificationMd();

		if (is2dClassification) {
			classes = new ClassInfo[ids.length];
			classesArray = new ArrayList<ClassInfo>();
			loadClassesInfo();
		}

		if (isRotSpectraMd() && mode == Mode.GALLERY_ROTSPECTRA) {
			if (zoom == 0)
				zoom = 100;
			return;
		}

		if (!md.isColumnFormat()) {
			mode = Mode.TABLE_MD;
			if (zoom == 0)
				zoom = 100;
		}

		if (isGalleryMode())
			mode = Mode.GALLERY_MD;

		if (hasRenderLabel()) {
			int renderLabel = ciFirstRender.getLabel();
			ImageGeneric image = null;
			String imageFn;
			// Try to find at least one image to render
			// and take dimensions from that
			for (int i = 0; i < ids.length && image == null; ++i) {
				imageFn = Filename.findImagePath(
						md.getValueString(renderLabel, ids[i]), filename, true);
				// DEBUG.printFormat("imageFn1: %s", imageFn);
				// imageFn = Filename.fixPath(md.getValueString(renderLabel,
				// ids[i]), filename, false);
				// DEBUG.printFormat("imageFn2: %s", imageFn);
				// if (imageFn != null){
				if (imageFn != null) {
					try {
						image = new ImageGeneric(imageFn);
					} catch (Exception e) {
						image = null;
					}
				}
				break;
			}
			if (image != null) { // Image file was found to render
				if (zoom == 0) { // if default value, calculate zoom
									// If in micrograph mode, reduce the
									// MAX_SIZE constant
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

					for (int i = 0; i < numberOfVols; ++i) {
						volumes[i] = md.getValueString(
								ciFirstRender.getLabel(), ids[i]);
					}
					commonVolPrefix = XmippStringUtils
							.commonPathPrefix(volumes);

					if (selectedVolFn.isEmpty())
						selectVolume(volumes[0]);

				}
				image.destroy();
			} else
				zoom = 100; // Render missing image icon at zoom 100
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
			int inputRenderLabel = MDLabel.MDL_UNDEFINED;

			if (!parameters.renderLabel.equalsIgnoreCase("first")) {
				inputRenderLabel = MetaData.str2Label(parameters.renderLabel);
			}

			for (int i = 0; i < lab.length; ++i) {
				ci = new ColumnInfo(lab[i]);
				if (labels != null) {
					for (ColumnInfo ci2 : labels)
						if (ci.label == ci2.label)
							ci.updateInfo(ci2);
				} else if (ci.allowRender)
					ci.render = globalRender;
				newLabels.add(ci);
				if (inputRenderLabel == lab[i] && ci.allowRender) {
					ciFirstRender = ci;
					if (ci.visible)
						ciFirstRenderVisible = ci;
				}
				if ((ciFirstRender == null || ci.getLabel() == MDLabel.MDL_IMAGE)
						&& ci.allowRender)// favor mdl_image over mdl_micrograph
				{
					ciFirstRender = ci;
				}
				if ((ciFirstRenderVisible == null || ci.getLabel() == MDLabel.MDL_IMAGE)
						&& ci.allowRender && ci.visible)
					ciFirstRenderVisible = ci;
			}
			if (ciFirstRenderVisible != null) {
				ciFirstRender = ciFirstRenderVisible;
			}
			// Add MDL_ENABLED if not present
			if (!md.containsLabel(MDLabel.MDL_ENABLED)
					&& (md.containsLabel(MDLabel.MDL_IMAGE) || md
							.containsLabel(MDLabel.MDL_MICROGRAPH))) {
				newLabels.add(0, new ColumnInfo(MDLabel.MDL_ENABLED));
				md.addLabel(MDLabel.MDL_ENABLED);
				for (long id : ids)
					md.setEnabled(true, id);
				// hasMdChanges = true;
			}

			labels = newLabels;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}// function loadLabels

	/** Read metadata and store ids */
	private void readMetadata(String fn) {
		try {
			hasMdChanges = false;
			hasClassesChanges = false;
			md.read(fn);

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			md = null;
			ids = null;
		}
	}

	/**
	 * Sort the metadata by a given column. The sort could be ascending or
	 * descending
	 */
	public void sortMd(int col, boolean ascending) {
		try {
			md.sort(getLabelFromCol(col), ascending);
			clearSelection();
			hasMdChanges = true;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	void clearSelection() {
		for (int i = 0; i < selection.length; ++i)
			selection[i] = false;

	}

	/** Reload current metadata from file */
	public void readMd() {
		if (filename != null)
			readMetadata(getMdFilename());
	}

	/** Select one of the blocks */
	public void selectBlock(String block) {
		selectedBlock = block;
		selectedVolFn = ""; // Set to empty string to get the first vol
		readMd();
	}

	public ImageGalleryTableModel createModel() {
		try {
			switch (mode) {
			case GALLERY_VOL:
				return new VolumeGalleryTableModel(this);
			case GALLERY_MD:
				if (md.size() > 0 && hasRenderLabel())
					return new MetadataGalleryTableModel(this);
				// else fall in the next case
			case TABLE_MD:
				mode = Mode.TABLE_MD; // this is necessary when coming from
				// previous case
				if (!md.isColumnFormat())
					return new MetadataRowTableModel(this);
				if (md.containsMicrographsInfo())
					return new MicrographsTableModel(this);
				return new MetadataTableModel(this);
			case GALLERY_ROTSPECTRA:
				return new RotSpectraGalleryTableModel(this);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public int getNumberOfBlocks() {
		return mdBlocks != null ? mdBlocks.length : 0;
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
	public boolean allowGallery() {
		return hasRenderLabel() || isRotSpectraMd();
	}

	// some mode shortcuts
	public boolean isGalleryMode() {
		return mode == Mode.GALLERY_MD || mode == Mode.GALLERY_VOL
				|| mode == Mode.GALLERY_ROTSPECTRA;
	}

	public boolean isVolumeMode() {
		return mode == Mode.GALLERY_VOL;
	}

	public boolean isTableMode() {
		return mode == Mode.TABLE_MD;
	}

	/** Return true if the underlying metadata is in row format */
	public boolean isColumnFormat() {
		return md.isColumnFormat();
	}

	public boolean isRotSpectraMode() {
		return mode == Mode.GALLERY_ROTSPECTRA;
	}

	public boolean isMicrographsMode() {
		return md.containsMicrographsInfo();
	}

	// utility function to change of mode
	public void changeMode() {
		if (isGalleryMode()) {
			mode = Mode.TABLE_MD;
			if (selection.length < ids.length) //This can happen when in volume mode, that changes the selection array
				selection = new boolean[ids.length];
		}
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
		selectedVolFn = vol; // FIXME: Check it is valid
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
			if (isVolumeMode() || !md.containsLabel(MDLabel.MDL_ENABLED))
				return true;
			return md.getEnabled(ids[index]);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return true;
	}

	/** Set enabled state */
	public void setEnabled(int index, boolean value) {
		try {
			if (!isVolumeMode()) { // slices in a volume are always enabled
				md.setEnabled(value, ids[index]);
				hasMdChanges = true;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/** Set all values coming from a row md */
	public void setRow(MDRow mdRow, long objId) {
		md.setRow(mdRow, objId);
		setMdChanges(true);
	}

	/** This is only needed for metadata table galleries */
	public boolean isFile(ColumnInfo ci) {
		try {
			return MetaData.isPathField(ci.getLabel());
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}

	public boolean isFile(int col) {
		return isFile(labels.get(col));
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

	public boolean is2DClassificationMd() {
		return is2dClassification;
	}

	/** Return true if current metadata comes from 2d classification */
	public boolean checkifIs2DClassificationMd() {
		try {
			if (!selectedBlock.startsWith("classes")
					|| !(md.containsLabel(MDLabel.MDL_REF) && md
							.containsLabel(MDLabel.MDL_CLASS_COUNT)))
				return false;
			for (long id : ids) {
				int ref = md.getValueInt(MDLabel.MDL_REF, id);
				long count = md.getValueLong(MDLabel.MDL_CLASS_COUNT, id);
				String s = Filename.getClassBlockName(ref);
				if (count > 0 && !containsBlock(s)) {
					// DEBUG.printFormat("2Dclass: for ref: %d, no block '%s'",
					// ref, s);
					return false;
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return true;
	}

	/** Get the assigned class of some element */
	public ClassInfo getItemClassInfo(int index) {
		if (is2dClassification && index < classes.length) {
			return classes[index];
		}
		return null;
	}

	/** Set item class info in md */
	private void setItemClassInfo(long id, ClassInfo cli) {
		String comment = "None";
		int ref2 = -1;
		int color = -1;
		try {
			if (cli != null) {
				ref2 = cli.index + 1;
				color = cli.getColor().getRGB();
				comment = cli.getComment();
			}
			md.setValueInt(MDLabel.MDL_REF2, ref2, id);
			md.setValueString(MDLabel.MDL_KEYWORDS, comment, id);
			md.setValueInt(MDLabel.MDL_COLOR, color, id);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	/** Set the class of an element */
	public void setItemClass(int index, ClassInfo cli) {
		hasClassesChanges = true;
		classes[index] = cli;
		long id = ids[index];
		setItemClassInfo(id, cli);
	}

	public ClassInfo getClassInfo(int classNumber) {
		return classesArray.get(classNumber);
	}

	/**
	 * Compute and update the number of classes and images assigned to this
	 * superclass
	 */
	public void updateClassesInfo() {
		try {
			int i = 0;
			for (ClassInfo cli : classesArray) {
				cli.numberOfClasses = 0;
				cli.numberOfImages = 0;
				cli.index = i++;
			}
			i = 0;
			for (long id : ids) { // iterate over all references
				long count = md.getValueLong(MDLabel.MDL_CLASS_COUNT, id);

				ClassInfo cli = getItemClassInfo(i);
				if (cli != null) {
					cli.numberOfClasses += 1;
					cli.numberOfImages += count;
					hasMdChanges = true;
				}
				setItemClassInfo(id, cli);
				++i;
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}// function upateClassesInfo

	/** Load classes structure if previously stored */
	public void loadClassesInfo() {
		try {
			if (md.containsLabel(MDLabel.MDL_REF2)) {
				long id;
				int ref2;
				ClassInfo cli;

				for (int i = 0; i < ids.length; ++i) {
					id = ids[i];
					ref2 = md.getValueInt(MDLabel.MDL_REF2, id);

					if (ref2 > 0) {
						cli = null;

						for (ClassInfo cli2 : classesArray)
							if (cli2.index == ref2) {
								cli = cli2;
								break;
							}

						if (cli == null) {
							String comment = md.getValueString(
									MDLabel.MDL_KEYWORDS, id);
							int color = md.getValueInt(MDLabel.MDL_COLOR, id);
							cli = new ClassInfo(comment, new Color(color));
							cli.index = ref2;
							classesArray.add(cli);
						}
						classes[i] = cli;
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}// function loadClassesInfo

	/** Return the number of selected elements */
	public int getSelectionCount() {
		int count = 0;
		if (!isVolumeMode()) {
			for (int i = 0; i < ids.length; ++i)
				if (selection[i])
					++count;
		}
		return count;
	}

	/** Create a metadata just with selected items */
	public MetaData getSelectionMd() {
		MetaData selectionMd = null;
		if (!isVolumeMode()) {
			long[] selectedIds = new long[getSelectionCount()];
			int count = 0;
			for (int i = 0; i < ids.length; ++i)
				if (selection[i])
					selectedIds[count++] = ids[i];
			try {
				selectionMd = new MetaData();
				selectionMd.importObjects(md, selectedIds);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return selectionMd;
	}

	/**
	 * Compute the metadatas
	 */
	public MetaData[] getClassesMd() {
		try {
			if (!classesArray.isEmpty()) {
				updateClassesInfo();
				// Count the number of non-empty classes
				MetaData[] mds = new MetaData[classesArray.size() + 1];
				mds[0] = new MetaData();
				MetaData mdAux = mds[0];
				int i = 0;
				long id;
				// Md for classes block
				for (ClassInfo cli : classesArray) {
					id = mdAux.addObject();
					mdAux.setValueInt(MDLabel.MDL_REF, ++i, id);
					mdAux.setValueLong(MDLabel.MDL_CLASS_COUNT,
							cli.numberOfImages, id);
					mdAux.setValueString(MDLabel.MDL_KEYWORDS,
							cli.getComment(), id);
					mds[i] = new MetaData();
				}
				i = 0;
				// Fill the classX_images blocks
				for (i = 0; i < ids.length; ++i) {
					ClassInfo cli = getItemClassInfo(i);
					if (cli != null) {
						id = ids[i];
						md.setValueInt(MDLabel.MDL_REF2, cli.index + 1, id);
						md.setValueString(MDLabel.MDL_KEYWORDS,
								cli.getComment(), id);
						mdAux = getClassImages(i);
						if (mdAux != null)
							mds[cli.index + 1].unionAll(mdAux);
					}
				}
				return mds;
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	/** Get the metadata with assigned images to this classes */
	public MetaData getClassImages(int index) {
		try {
			long id = ids[index];
			int ref = md.getValueInt(MDLabel.MDL_REF, id);
			String blockName = Filename.getClassBlockName(ref);
			if (containsBlock(blockName)) {
				return new MetaData(blockName + Filename.SEPARATOR + filename);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/** Return true if current metadata is a rotspectra classes */
	public boolean isRotSpectraMd() {
		if (filename != null) {
			String fnVectors = filename.replace("classes", "vectors");
			String fnVectorsData = fnVectors.replace(".xmd", ".vec");
			if (is2DClassificationMd() && Filename.exists(fnVectors)
					&& Filename.exists(fnVectorsData))
				return true;
		}
		return false;
	}

	/** Check if a block is present, ignore case */
	public boolean containsBlock(String block) {
		if (mdBlocks != null)
			for (String b : mdBlocks)
				if (b.equalsIgnoreCase(block))
					return true;
		return false;
	}

	/**
	 * Take an index counting only visible columns and translate into the
	 * general column index
	 * 
	 * @param col
	 *            column index in visible counting
	 * @return column index in general counting
	 */
	public int getVisibleColumnIndex(int col) {
		int visibleIndex = 0;
		for (int i = 0; i < labels.size(); i++)
			if (labels.get(i).visible) {
				if (col == visibleIndex)
					return i;
				visibleIndex++;
			}
		return -1;
	}

	public int getLabelFromCol(int col) {
		return labels.get(col).getLabel();
	}

	public ColumnInfo getColumnInfo(int col) {
		return labels.get(col);
	}

	public String getValueFromCol(int index, int col) {
		if (!isColumnFormat()) {
			col = index;
			index = 0;
		}
		return getValueFromCol(index, labels.get(col));
	}

	public String getValueFromCol(int index, ColumnInfo ci) {
		try {
			return md.getValueString(ci.getLabel(), ids[index]);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public String getValueFromLabel(int index, int label) {
		try {
			return md.getValueString(label, ids[index]);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public void setValueToCol(int index, ColumnInfo ci, String value) {
		try {
			md.setValueString(ci.getLabel(), value, ids[index]);
			setMdChanges(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/** Delete from metadata selected items */
	public void removeSelection() throws Exception {
		for (int i = 0; i < ids.length; ++i) {
			if (selection[i]) {
				md.removeObject(ids[i]);
				hasMdChanges = true;
			}
		}
	}

	/** Add a new class */
	public void addClass(ClassInfo ci) {
		classesArray.add(ci);
		hasClassesChanges = true;
	}

	/** Remove a class from the selection */
	public void removeClass(int classNumber) {
		ClassInfo cli = getClassInfo(classNumber);
		for (int i = 0; i < ids.length; ++i)
			if (getItemClassInfo(i) == cli)
				setItemClass(i, null);
		classesArray.remove(classNumber);
		hasClassesChanges = true;
	}

	public boolean hasMdChanges() {
		return hasMdChanges;
	}

	public void setMdChanges(boolean value) {
		hasMdChanges = value;
	}

	public boolean hasClassesChanges() {
		return hasClassesChanges;
	}

	public boolean hasMicrographParticles() {
		return md.containsMicrographParticles();
	}

	public String getFileName() {
		return filename;

	}

	public String getBlock(int index) {
		int size = getNumberOfBlocks();
		if (size > 0 && index >= 0 && index < size)
			return mdBlocks[index];
		return null;
	}

	public MetaData getImagesMd() {
		int idlabel = getRenderLabel();
		if (md == null)
			return null;
		if (!md.containsLabel(idlabel))
			return null;

		MDRow mdRow = new MDRow();
		MetaData imagesmd = new MetaData();
		int index = 0;
		String imagepath;
		long id2;
		// md.print();
		for (long id : md.findObjects()) {
			if (isEnabled(index)) {
				imagepath = md.getValueString(idlabel, id, true);
				if (imagepath != null && ImageGeneric.exists(imagepath)) {
					id2 = imagesmd.addObject();
					if (useGeo) {
						md.getRow(mdRow, id);
						mdRow.setValueString(idlabel, imagepath);
						imagesmd.setRow(mdRow, id2);
					} else
						imagesmd.setValueString(idlabel, imagepath, id2);
				}
			}
			index++;
		}
		mdRow.destroy();
		return imagesmd;
	}

	public String getFileInfo() {
		File file = new File(getFileName());

		String fileInfo = "Path: " + file.getAbsolutePath() + "\n\n";

		fileInfo += "File Name: " + file.getName() + "\n" + "Last Modified: "
				+ new Date(file.lastModified()) + "\n"

				+ "Size: " + Filename.humanReadableByteCount(file.length());
		return fileInfo;
	}
}// class GalleryDaa
