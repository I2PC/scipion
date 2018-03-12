/**
 * *************************************************************************
 * Authors: J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * All comments concerning this program package may be sent to the e-mail
 * address 'xmipp@cnb.csic.es'
 **************************************************************************
 */
package xmipp.viewer.models;

import ij.ImagePlus;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import xmipp.ij.commons.Geometry;
import xmipp.ij.commons.Tool;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.EllipseCTF;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MDRow;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.Params;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippStringUtils;
import xmipp.viewer.ctf.CTFAnalyzerJFrame;
import xmipp.viewer.ctf.CTFRecalculateImageWindow;
import xmipp.viewer.ctf.EstimateFromCTFTask;
import xmipp.viewer.ctf.TasksEngine;
import xmipp.viewer.scipion.ScipionGalleryData;
import xmipp.viewer.scipion.ScipionMetaData;
import xmipp.viewer.windows.AddObjectJDialog;
import xmipp.viewer.windows.GalleryJFrame;
import xmipp.viewer.windows.SaveJDialog;

/**
 * This class will serve to store important data about the gallery
 */
public class GalleryData {


    protected HashMap<String, ColumnInfo> displaycis;
    protected MetaData md;
    protected long[] ids;
    protected String[] mdBlocks = null;
    protected String selectedBlock;
    // The following is only used in VolumeGallery mode
    protected String selectedVolFn = "";
    //protected String commonVolPrefix = "";

    protected List<ColumnInfo> labels = null;
    // First label that can be rendered
    protected ColumnInfo ciFirstRender = null;
    protected Integer zoom;
    protected String filename;
    protected int resliceView;
    protected Mode mode;
    protected boolean showLabel = false;
    protected boolean renderImages;
    public final Params parameters;


    // flag to perform global normalization
    protected boolean normalize = false;
    // flag to use geometry info
    protected boolean useGeo;
    // flag to wrapping
    protected boolean wrap;
    // flag to check if is 2d classification
    protected boolean isClassification = false;
    protected int refLabel;
    // Store the selection state for each item
    

    // Array with all ClassInfo
    protected ArrayList<ClassInfo> classesArray;
    // ClassInfo reference for each element
    protected ClassInfo[] classes;
    // Flags to check if md or classes has changed
    protected boolean hasMdChanges, hasClassesChanges;
    protected GalleryJFrame window;
    
    protected String[] displayLabels;
    protected String[] sortby;
    
    protected boolean isVolumeMd;
    protected Integer rows, columns;

    
    
    protected String[] renderLabels;
    protected String renderLabel = "first";
    protected String[] visibleLabels;
    protected String[] orderLabels;
    protected boolean inverty;

    public boolean isObjectCmd(String cmd) {
        for(String objCmd: parameters.objectCommands)
            if(objCmd.equals(cmd))
                return true;
        return false;
    }

    public boolean isVolumeMd()
    {
    	return isVolumeMd;
    }
    
    public Integer getModelRows()
    {
        return rows;
    }
    
    public Integer getModelColumns()
    {
        return columns;
    }
    
    public void runObjectCommand(int index, String objectCommand) {
        try {
            String[] cmd = new String[]{objectCommand};
            System.out.println(Arrays.toString(cmd));
        } catch (Exception ex) {
            Logger.getLogger(GalleryData.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public String getTmpFile(String suffix) {
        String ext = XmippStringUtils.getFileExtension(filename);
        String ext2 = suffix + ext;
        String tmpfile = filename;
        if(!filename.endsWith(ext2))
            tmpfile = filename.replace(ext, ext2);
        return tmpfile;
    }

    public void selectVolumeAt(int selectedIndex) {
    	String vol = getVolumeAt(selectedIndex);
    	String file = vol;
    	if (vol.contains(":"))
    		file = vol.substring(0, vol.lastIndexOf(":"));
    	if (vol.contains("@"))
    		file = file.substring(file.lastIndexOf("@") + 1);
    	if(new File(file).exists())
    	{
    		selectedVolFn = vol;
    	}
    	else
    		throw new IllegalArgumentException(XmippMessage.getPathNotExistsMsg(file));
    }

    public boolean isAutoAdjust()
    {
        return rows == null && columns == null;
    }

    public void setModelDim(Integer rows, Integer cols) {
    	if(cols != null && cols == 0)
    		throw new IllegalArgumentException(XmippMessage.getIllegalValueMsg("columns", 0));
        this.rows = rows;
        this.columns = cols;
    }

    public boolean isScipionInstance() {
        return this instanceof ScipionGalleryData;
    }

    public boolean isInvertY() {
        return inverty;
    }

    public String getChimeraProjectionCmd(int row) {
        String rot = getValueFromLabel(row, MDLabel.MDL_ANGLE_ROT);
        String tilt = getValueFromLabel(row, MDLabel.MDL_ANGLE_TILT);
        String psi = getValueFromLabel(row, MDLabel.MDL_ANGLE_PSI);
        String command = String.format("rotate %s %s %s", rot, tilt, psi);
        return command;
    }


    public enum Mode {

        GALLERY_MD, GALLERY_VOL, TABLE_MD
    };

    // define min and max render dimensions
    public static int MIN_SIZE = 16;
    public static int MAX_SIZE = 256;

    public GalleryData(GalleryJFrame window, String fn, Params parameters) 
    {
        this(window, parameters, new MetaData(fn));
    }
	// max dimension allowed to render images
    /**
     * The constructor receive the filename of a metadata The metadata can also
     * be passed, if null, it will be readed from filename
     *
     * @param jFrameGallery
     */
    public GalleryData(GalleryJFrame window, Params parameters, MetaData md) {
        this.window = window;
        this.parameters = parameters;
        
        try {
            selectedBlock = "";
            this.md = md;
            setFileName(md.getPath());
            zoom = parameters.zoom;
            resliceView = parameters.resliceView;
            if(parameters.rows > 0)
                rows = parameters.rows;
            if(parameters.columns > 0)
                columns = parameters.columns;
            loadMd();

        } catch (Exception e) {
            e.printStackTrace();
            md = null;
        }

    }// constructor GalleryData
    
    protected void readMdParameters()
    {
        renderLabels = null;
        renderLabel = "first";
        visibleLabels = null;
        orderLabels = null;
        sortby = null;
        
        useGeo = wrap = false;
        displayLabels = null;

        if (mode == null)
            mode = Mode.GALLERY_MD;
        this.renderImages = true;
        displaycis = new HashMap<String, ColumnInfo>();
        this.inverty = parameters.inverty;

        if(parameters.getBlock() == null)
            parameters.setBlock(selectedBlock);//Identifies parameters with first block loaded

        if (parameters.mode.equalsIgnoreCase(Params.OPENING_MODE_METADATA))
            mode = Mode.TABLE_MD;
        else if (parameters.mode.equalsIgnoreCase(Params.OPENING_MODE_GALLERY))
            mode = Mode.GALLERY_MD;

        if(parameters.getBlock().equals(selectedBlock))
        {
            setRenderLabels(parameters.renderLabels);
            setRenderLabel(parameters.getRenderLabel());
            setVisibleLabels(parameters.visibleLabels);
            setOrderLabels(parameters.orderLabels);
            useGeo = parameters.useGeo;
            wrap = parameters.wrap;
            displayLabels = parameters.getDisplayLabels();
            this.renderImages = parameters.renderImages;
            
        }
    }

    public List<ColumnInfo> getColumns() {
        return labels;
    }

    public void setRenderColumn(String key) {
            if(key.equalsIgnoreCase("none"))
                ciFirstRender = null;
            for(ColumnInfo ci: labels)
                if(ci.labelName.equals(key))
                    ciFirstRender = ci;
            renderLabel = ciFirstRender.labelName;
	}

    public ColumnInfo getRenderColumn() {
        return ciFirstRender;
    }

    /**
     * Return the name of the selected md block
     */
    public String getMdFilename() {
        if (selectedBlock.isEmpty()) {
            return filename;
        }
        return String.format("%s@%s", selectedBlock, filename);
    }// function getMdFilename

    public String getMdSaveFileName() {
        if (filename == null) {
            return null;
        }

        String savefn;
        if (selectedBlock.isEmpty()) {
            savefn = filename;
        } else {
            savefn = String.format("%s@%s", selectedBlock, filename);
        }
        String ext;
        if (savefn.contains(".")) {
            ext = savefn.substring(savefn.lastIndexOf("."));
            if (ext.equals(".stk")) {
                return savefn.replace(ext, ".xmd");
            }
        } else {
            savefn = savefn + ".xmd";
        }
        return savefn;
    }

    public void setFileName(String file) {
        filename = file;
        if (file != null) {
            if (Filename.hasPrefix(file)) {
                if (Filename.isMetadata(file)) {
                    selectedBlock = Filename.getPrefix(file); 
                    filename = Filename.getFilename(file);
                }
            }
            if (Filename.exists(filename)) {
                mdBlocks = MetaData.getBlocksInMetaDataFile(filename);
                if (mdBlocks.length >= 1 && selectedBlock.isEmpty()) {
                    selectedBlock = mdBlocks[0];
                }
            }
        }
    }

    public void setDisplayLabel(String key, boolean selected) {
        if(!selected)
        {
            displaycis.remove(key);
            return;
        }
        
        for(ColumnInfo ci: labels)
            if(ci.labelName.equals(key))
            {
                displaycis.put(key, ci);
                break;
            }
    }

    /**
     * Load contents from a metadata already read
     */
    public void loadMd() throws Exception {
        readMdParameters();
        ids = md.findObjects();
        loadLabels();
        
        if (!containsGeometryInfo()) {
            useGeo = false;
            wrap = false;
        }
        
        isClassification = checkifIsClassificationMd();

        if (isClassification) {
            classes = new ClassInfo[ids.length];
            classesArray = new ArrayList<ClassInfo>();
            loadClassesInfo();
        }

        
        if (!md.isColumnFormat() ) {
            mode = Mode.TABLE_MD;
            renderImages = false;
        }

        if (isGalleryMode()) {
            mode = Mode.GALLERY_MD;
        }
        isVolumeMd = false;
        if (hasRenderLabel()) {
            int renderLabel = ciFirstRender.label;
            ImageGeneric image = null;
            String imageFn = null;


			// Try to find at least one image to render
            // and take dimensions from that
            for (int i = 0; i < ids.length && image == null; ++i) {
            	if(i == 100)//after 100 items if there is no image available break 
            		break;
                imageFn = md.getValueString(renderLabel, ids[i]);
                if(imageFn != null)
                    imageFn = Filename.findImagePath(imageFn , filename, true);
                if (imageFn != null) {
                    try {
                        image = new ImageGeneric(imageFn);
                    } catch (Exception e) {
                        image = null;
                    }
                    break;
                }
                
            }
            if (image != null) { // Image file was found to render
                if (zoom == null) { // if default value, calculate zoom
                    // If in micrograph mode, reduce the
                    // MAX_SIZE constant
                    if (md.containsMicrographsInfo()) {
                        MAX_SIZE /= 2;
                    }
                    zoom = getDefaultZoom(image.getXDim());
                    
                }
                if (image.isVolume()) { // We are assuming all are volumes
                    // or images, don't mix it
                    if (isGalleryMode()) 
                        mode = Mode.GALLERY_VOL;
                    
                    isVolumeMd = true;
                    if (selectedVolFn.isEmpty()) {
                        selectedVolFn = imageFn;
                    }

                }
                image.destroy();
            } 
        } else {
            // force this mode when there aren't render label
            mode = Mode.TABLE_MD;
        }
        if(parameters.sortby != null && parameters.getBlock().equals(selectedBlock))
        {
            ColumnInfo sortci = getColumnInfo(parameters.sortby[0]);
            boolean asc = parameters.sortby.length == 1 || parameters.sortby[1].equalsIgnoreCase("asc");
            if(sortci != null)
                sortMd(sortci, asc);
        }
        
    }// function loadMd
    
 
    public static int getDefaultZoom(int xdim)
    {
        int x = Math.min(Math.max(xdim, MIN_SIZE), MAX_SIZE);
        float scale = (float) x / xdim;
        int zoom = (int) Math.ceil(scale * 100);
        return zoom;
    }

    public ColumnInfo getColumnInfo(String labelName) {
        for (ColumnInfo ci : labels) {
            if (ci.labelName.equals(labelName)) {
                return ci;
            }
        }
        return null;
    }
        
               
    public boolean isDisplayLabel()
    {
        return displaycis != null;
    }

   

    public String getDisplayLabel(long id)
    {
        String label = "";
        
        for(ColumnInfo ci: displaycis.values())
            label += ci.labelName + "=" + md.getValueString(ci.label, id) + ", ";
        if(!label.isEmpty())
            label = label.substring(0, label.length() - 2);
        return label;
    }

    /**
     * Load labels info in md, try to keep previous settings of render and
     * visible on same columns
     */
    public void loadLabels() {

        ColumnInfo ci;
        try {
            int[] labelids = md.getActiveLabels();
            ArrayList<ColumnInfo> newLabels = new ArrayList<ColumnInfo>(
                    labelids.length);
            ciFirstRender = null;
            ColumnInfo ciFirstRenderVisible = null;
            String inputRenderLabel = "";

            if (!renderLabel.equalsIgnoreCase("first")) {
                inputRenderLabel = renderLabel;
            }
            String sampleImage;
            ColumnInfo auxRender = null;
            for (int i = 0; i < labelids.length; ++i) {
                ci = initColumnInfo(labelids[i]);

                ci.render = isRenderLabel(ci);
                if(ci.render)
                    auxRender = ci;
                ci.visible = isVisibleLabel(ci);
                newLabels.add(ci);
                if (inputRenderLabel.equals(ci.labelName) && ci.render) {//render label specified and included on renders
                    ciFirstRender = ci;
                    if (ci.visible) {
                        ciFirstRenderVisible = ci;
                    }
                }
                

                if (ciFirstRender == null && ci.allowRender)
                {
                    sampleImage = getSampleImage(ci);
                    if(sampleImage != null)
                    {
                        ciFirstRender = ci;
                        if (ciFirstRenderVisible == null && ci.visible) 
                            ciFirstRenderVisible = ci;
                    }
                }
            }
            if(ciFirstRender == null)
                ciFirstRender = ciFirstRenderVisible = auxRender;//if there are no images found render none image on gallery mode
            if (ciFirstRenderVisible != null) {
                ciFirstRender = ciFirstRenderVisible;
            }
            
            // Add MDL_ENABLED if not present
            if (!md.containsLabel(MDLabel.MDL_ENABLED) && (md.containsLabel(MDLabel.MDL_IMAGE) || md.containsLabel(MDLabel.MDL_MICROGRAPH))) {
                newLabels.add(0, new ColumnInfo(MDLabel.MDL_ENABLED));
                md.addLabel(MDLabel.MDL_ENABLED);
                for (long id : ids) {
                    md.setEnabled(true, id);
                }
                // hasMdChanges = true;
            }
            

            labels = newLabels;
            orderLabels();
            if(displayLabels != null)
                for(String label: displayLabels)
                    setDisplayLabel(label, true);
//////            System.out.printf("render: %s %s \n", ciFirstRender, ciFirstRenderVisible);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }// function loadLabels
    
    public String getSampleImage(ColumnInfo ci)
    {
        String imageFn, mddir = md.getBaseDir();
        for (int i = 0; i < ids.length; ++i)
        {
        	if( i == 100)
        		break;
            imageFn = getValueFromLabel(i, ci.label);
            if(imageFn != null)
            {
                imageFn = Filename.findImagePath(imageFn, mddir, true);
                if (imageFn != null && Filename.exists(imageFn))
                       return imageFn;
            }
        }
        return null;
    }
    
    public ColumnInfo initColumnInfo(int label)
    {
        return new ColumnInfo(label);
    }

    /**
     * Read metadata and store ids
     */
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
    public void sortMd(ColumnInfo sortci, boolean asc) {
        try {
            if(sortby == null)
                sortby = new String[2];
            sortby[0] = sortci.labelName;
            sortby[1] = (asc)? "ASC": "DESC";
            md.sort(sortci.label, asc);
            hasMdChanges = true;
            ids = md.findObjects();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    

   

    /**
     * Reload current metadata from file
     */
    public void readMd() {
        if (filename != null) {
            readMetadata(getMdFilename());
        }
    }

    /**
     * Select one of the blocks
     */
    public void selectBlock(String block) {
        selectedBlock = block;
        selectedVolFn = ""; // Set to empty string to get the first vol
        readMd();
    }

    /**
     * defines how each row is rendered on ImageGalleryTableModel
     *
     * @return
     */
    public ImageGalleryTableModel createModel(boolean[] selection) {
        try {
            switch (mode) {
                case GALLERY_VOL:
                    return new VolumeGalleryTableModel(this);
                case GALLERY_MD:
                    if (md.size() > 0 && hasRenderLabel()) {
                    	
                        return new MetadataGalleryTableModel(this, selection);
                    }
                // else fall in the next case
                case TABLE_MD:
                    mode = Mode.TABLE_MD; // this is necessary when coming from
                    // previous case
                    if (!md.isColumnFormat()) {
                        return  new MetadataRowTableModel(this);
                    }

                    return new MetadataTableModel(this, selection);
                
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    public int getNumberOfBlocks() {
        return mdBlocks != null ? mdBlocks.length : 0;
    }

   

    /**
     * Return the mode of the gallery
     */
    public Mode getMode() {
        return mode;
    }

    /**
     * Return true if there is a renderizable label in the metadata
     */
    public boolean hasRenderLabel() {
        return ciFirstRender != null;
    }

    /**
     * Return the label that is used for rendering
     */
    public int getRenderLabel() {
        return ciFirstRender.label;
    }

    /**
     * Return true if the gallery mode is allowed
     */
    public boolean allowGallery() {
        return hasRenderLabel();
    }

    // some mode shortcuts
    public boolean isGalleryMode() {
        return mode == Mode.GALLERY_MD || mode == Mode.GALLERY_VOL;
                
    }

    public boolean isVolumeMode() {
        return mode == Mode.GALLERY_VOL;
    }

    public boolean isTableMode() {
        return mode == Mode.TABLE_MD;
    }

    /**
     * Return true if the underlying metadata is in row format
     */
    public boolean isColumnFormat() {
        return md.isColumnFormat();
    }
    

    public boolean isMicrographsMode() {
        return md.containsMicrographsInfo();
    }

    // utility function to change of mode
    public void changeMode() {
        if (isGalleryMode()) {
            mode = Mode.TABLE_MD;
           
        } 
        else if ( isVolumeMd) 
            mode = Mode.GALLERY_VOL;
        else
            mode = Mode.GALLERY_MD;
    }

    /**
     * following function only should be used in VolumeGallery mode
     */
    public String getVolumeAt(int index) {
            return md.getValueString(ciFirstRender.label, ids[index]);
    }

    public boolean containsGeometryInfo() 
    {
        return containsGeometryInfo("2D");
    }

    // Check if the underlying data has geometrical information
    public boolean containsGeometryInfo(String type) 
    {
        try {
            return md.containsGeometryInfo(type);
        } catch (Exception e) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Check if an item is enabled or not
     */
    public boolean isEnabled(int index) {
        try {
            if (isVolumeMode()) {
                return true;
            }
            return md.getEnabled(ids[index]);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return true;
    }

    /**
     * Set enabled state
     */
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

    /**
     * Set all values coming from a row md
     */
    public void setRow(MDRow mdRow, long objId) {
        md.setRow(mdRow, objId);
        setMdChanges(true);
    }

    /**
     * This is only needed for metadata table galleries
     */
    public boolean isFile(ColumnInfo ci) {
        try {
            return MetaData.isPathField(ci.label);
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
        		
            return ci.render;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return false;
    }

    public boolean isClassificationMd() {
        return isClassification;
    }

    /**
     * Return true if current metadata comes from 2d classification
     */
    public boolean checkifIsClassificationMd() {
        try {
            boolean valid = selectedBlock.startsWith("classes")
                    && (md.containsLabel(MDLabel.MDL_REF) || md.containsLabel(MDLabel.MDL_REF3D))
                    && md.containsLabel(MDLabel.MDL_CLASS_COUNT);

            if (!valid) {
                return false;
            }

            refLabel = md.containsLabel(MDLabel.MDL_REF) ? MDLabel.MDL_REF : MDLabel.MDL_REF3D;

            for (long id : ids) {
                int ref = md.getValueInt(refLabel, id);
                long count = md.getValueLong(MDLabel.MDL_CLASS_COUNT, id);
                String s = Filename.getClassBlockName(ref);
                if (count > 0 && !containsBlock(s)) {
                    DEBUG.printFormat("2Dclass: for ref: %d, no block '%s'\n", ref, s);
                    return false;
                }
            }
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return true;
    }

    /**
     * Get the assigned class of some element
     */
    public ClassInfo getItemClassInfo(int index) {
        if (isClassification && index < classes.length) {
            return classes[index];
        }
        return null;
    }

    /**
     * Set item class info in md
     */
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

    /**
     * Set the class of an element
     */
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

    /**
     * Load classes structure if previously stored
     */
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

                        for (ClassInfo cli2 : classesArray) {
                            if (cli2.index == ref2) {
                                cli = cli2;
                                break;
                            }
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

    

    /**
     * Create a metadata just with selected items
     */
    public MetaData getSelectionMd(boolean[] selection) {
        MetaData selectionMd = null;
        if (!isVolumeMode()) {
            long[] selectedIds = new long[selection.length];
            int count = 0;
            for (int i = 0; i < selection.length; i++) {
                if (selection[i]) {
                    selectedIds[count++] = ids[i];
                }
            }
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
                        if (mdAux != null) {
                            mds[cli.index + 1].unionAll(mdAux);
                        }
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

    /**
     * Get the metadata with assigned images to this classes
     */
    public MetaData getClassImages(int index) {
        try {
            long id = ids[index];
            int ref = md.getValueInt(refLabel, id);
            String blockName = Filename.getClassBlockName(ref);
            if (containsBlock(blockName)) {
                return new MetaData(blockName + Filename.SEPARATOR + filename);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * Get all the images assigned to all selected classes
     */
    public MetaData getClassesImages(boolean[] selection) {
        MetaData mdImages = new MetaData();
        MetaData md;
        for (int i = 0; i < selection.length; i++) {
            if (selection[i]) {
                md = getClassImages(i);
                if (md != null) {
                    mdImages.unionAll(md);
                    md.destroy();
                }
            }
        }
        return mdImages;
    }

    public MetaData getEnabledClassesImages() {
        MetaData mdImages = new MetaData();
        MetaData md;
        for (int i = 0; i < ids.length; ++i) {
            if (isEnabled(i)) {
                md = getClassImages(i);
                if (md != null) {
                    mdImages.unionAll(md);
                    md.destroy();
                }
            }
        }
        return mdImages;
    }
    

    /**
     * Check if a block is present, ignore case
     */
    public boolean containsBlock(String block) {
        if (mdBlocks != null) {
            for (String b : mdBlocks) {
                if (b.equalsIgnoreCase(block)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Take an index counting only visible columns and translate into the
     * general column index
     *
     * @param col column index in visible counting
     * @return column index in general counting
     */
    public int getVisibleColumnIndex(int col) {
        int visibleIndex = 0;
        for (int i = 0; i < labels.size(); i++) {
            if (labels.get(i).visible) {
                if (col == visibleIndex) {
                    return i;
                }
                visibleIndex++;
            }
        }
        return -1;
    }

    public int getLabelFromCol(int col) {
        return labels.get(col).label;
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
            return md.getValueString(ci.label, ids[index]);
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
            md.setValueString(ci.label, value, ids[index]);
            setMdChanges(true);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Delete from metadata selected items
     */
    public void removeSelection(boolean[] selection) throws Exception {
        for (int i = 0; i < selection.length; i++) {
            if (selection[i]) {
                md.removeObject(ids[i]);
                hasMdChanges = true;
            }
        }
    }

    /**
     * Add a new class
     */
    public void addClass(ClassInfo ci) {
        classesArray.add(ci);
        hasClassesChanges = true;
    }

    /**
     * Remove a class from the selection
     */
    public void removeClass(int classNumber) {
        ClassInfo cli = getClassInfo(classNumber);
        for (int i = 0; i < ids.length; ++i) {
            if (getItemClassInfo(i) == cli) {
                setItemClass(i, null);
            }
        }
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
        if (size > 0 && index >= 0 && index < size) {
            return mdBlocks[index];
        }
        return null;
    }

    
    
    public MDRow[] getImages(MetaData md) {
        int idlabel = getRenderLabel();
        if (md == null) {
            return null;
        }
        if (!md.containsLabel(idlabel)) {
            return null;
        }

        MDRow mdRow = null;
        ArrayList<MDRow> imagesmd = new ArrayList<MDRow>();
        int index = 0;
        String imagepath;

        for (long id : md.findObjects()) {
            if (isEnabled(index)) {
                
                imagepath = md.getValueString(idlabel, id, true);
                if (imagepath != null && ImageGeneric.exists(imagepath)) {
                    mdRow = new MDRow();
                    
                    if (useGeo) {
                        
                        md.getRow(mdRow, id);
                        mdRow.setValueString(idlabel, imagepath);
                        
                    } else {
                        mdRow.setValueString(MDLabel.MDL_IMAGE, imagepath);
                    }
                    
                    imagesmd.add(mdRow);
                }
            }
            index++;
        }

        return imagesmd.toArray(new MDRow[]{});
    }

    public String getFileInfo() {
        File file = new File(getFileName());

        String fileInfo = "Path: " + file.getAbsolutePath() + "\n";

        fileInfo += "Last Modified: "
                + new Date(file.lastModified()) + "\n"
                + "Size: " + Filename.humanReadableByteCount(file.length());
        return fileInfo;
    }

    public void saveClassSelection(boolean[]selection, String path) {
        try {
            saveSelection(selection, "classes" + Filename.SEPARATOR + path, true);
            MetaData imagesmd;
            // Fill the classX_images blocks
                for (int i = 0; i < selection.length; i++) {
                    if (selection[i]) {
                        long id = ids[i];
                        int ref = md.getValueInt(MDLabel.MDL_REF, id);
                        String blockName = Filename.getClassBlockName(ref);
                        if (containsBlock(blockName)) {
                            imagesmd = new MetaData(blockName + Filename.SEPARATOR + filename);
                            imagesmd.writeBlock(blockName + Filename.SEPARATOR + path);
                            imagesmd.destroy();
                        }
                    }
                }

        } catch (Exception ex) {
            Logger.getLogger(GalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public boolean isCTFMd() {
        return md.isCTFMd();
    }

    /**
     * Save selected items as a metadata
     */
    public void saveSelection(boolean[] selection) throws Exception {

        SaveJDialog dlg = new SaveJDialog(window, "selection" + getFileExtension(), true);
        boolean save = dlg.showDialog();
        if (save) {
            boolean overwrite = dlg.isOverwrite();
            String path = dlg.getMdFilename();
            saveSelection(selection, path, overwrite);
        }

    }

    /**
     * Save selected items as a metadata
     */
    public void saveSelection(boolean[] selection, String path, boolean overwrite) throws Exception {
        MetaData md = getSelectionMd(selection);

        String file = path.substring(path.lastIndexOf("@") + 1, path.length());
        if (!new File(file).exists())// overwrite or append, save selection
        {
            md.write(path);
        } else {
            if (overwrite) {
                md.write(path);// overwrite with active block only, other
            } // blocks were dismissed
            else {
                md.writeBlock(path);// append selection
            }
        }

        md.destroy();
    }

    public int size() {
        return ids.length;
    }

    public Integer getZoom() {
        return zoom;
    }

    public boolean useGeo() {
        return useGeo;
    }

    public void fillConstant(int label, String value) {
        md.fillConstant(label, value);
    }

    public void fillLinear(int label, Double start, Double step) {
        md.fillLinear(label, start, step);
    }

    public void fillRandom(int label, String mode, double op1, double op2) {
        md.fillRandom(label, mode, op1, op2);
    }

    public void removeDisabled() {
        md.removeDisabled();
    }

    public void removeLabel(int label) {
        md.removeLabel(label);
    }

    public long[] getIds() {
        return ids;
    }

    public int[] getLabels() {
        return md.getActiveLabels();
    }

    public String getValueString(int label, long id) {
        return md.getValueString(label, id);
    }

    public boolean setValueString(int label, String newValue, long l) {
        return md.setValueString(label, newValue, l);
    }

    public Object getSelectedBlock() {
        return selectedBlock;
    }

    

    public String getSelVolumeFile() {
        return selectedVolFn;
    }

    public void setResliceView(int view) {
        resliceView = view;
    }

    public boolean isWrap() {
        return wrap;
    }

    public boolean renderImages() {
        return renderImages;
    }

    public int getResliceView() {
        return resliceView;
    }

    public boolean addObject() {
        AddObjectJDialog dlg = new AddObjectJDialog(window);
        if (dlg.showDialog()) {
            md.unionAll(dlg.md);
            return true;
        }
        return false;
    }

    public String[] getBlocks() {
        return mdBlocks;
    }

    public long getId(int i) {
        return ids[i];
    }

    public String createSortFile(String psdFile, int row) {

        MetaData mdRow = new MetaData();
        MDRow row2 = new MDRow();
        md.getRow(row2, getId(row));
        mdRow.setRow(row2, mdRow.addObject());
        String sortFn = psdFile.replace(".psd", ".xmd");
        mdRow.write(sortFn);
        mdRow.destroy();
        return sortFn;

    }

    public ArrayList<ClassInfo> getClasses() {
        return classesArray;
    }

    public List<ColumnInfo> getLabelsInfo() {
        return labels;
    }

    public MetaData getMetaDataRow() {
        return md.getMetaDataRow();
    }

    public String getLabel(long objId, int label) {
        try {
            if (isClassification) {
                int ref = md.getValueInt(MDLabel.MDL_REF, objId);
                long count = md.getValueLong(MDLabel.MDL_CLASS_COUNT, objId);
                return String.format("class %d (%d images)", ref, count);
            } else {
                return md.getValueString(label, objId);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    public String getFileExtension() {
        if (getFileName() == null) {
            return "";
        }
        return XmippStringUtils.getFileExtension(filename);
    }

    public void saveAll(String path, boolean overwrite) throws Exception {
        String from = getFileName();
        String blockto = path;
        String to;
        String block;
        if (blockto.contains("@")) {
            int sep = blockto.lastIndexOf("@");
            block = blockto.substring(0, sep);
            to = blockto.substring(sep + 1, blockto.length());
        } else {
            to = blockto;
            blockto = selectedBlock + "@" + blockto;
            block = selectedBlock;
        }

        if (from != null) {
            MetaData md;
            Hashtable<String, MetaData> mds = new Hashtable<String, MetaData>();
            for (String blockit : getBlocks()) {
                mds.put(blockit, getMetaData(blockit));
            }
            File file = new File(to);
            if (overwrite) {
                file.delete();
            }
            if (!file.exists() && file.getParentFile() != null) {
                file.getParentFile().mkdirs();
            }
            for (String blockit : getBlocks()) {
                md = mds.get(blockit);
                if (blockit.equals(selectedBlock)) {
                    saveMd(blockto, true, overwrite);
                } else {
                    md.writeBlock(blockit + "@" + to);
                }
                md.destroy();
            }
        } else {
            saveMd(blockto, true, overwrite);
        }

        setMdChanges(false);
        setFileName(to);
        if (blockto.contains("@")) {
            selectBlock(block);
        }
    }

    public void saveMd(String path, boolean saveall, boolean isoverwrite) throws Exception {
        try {

            if (path == null) {
                throw new IllegalArgumentException();
            }

            boolean overwritewithblock;
            String file;
            if (path.contains("@")) {
                file = path.substring(path.lastIndexOf("@") + 1);
            } else {
                file = path;
                path = selectedBlock + "@" + file;
            }

            File iofile = new File(file);
            if (!iofile.exists())// overwrite or append, save active
            {
                if (iofile.getParentFile() != null) {
                    iofile.getParentFile().mkdirs();
                }
                md.write(path);
            } else {

                overwritewithblock = isoverwrite && !saveall;
                if (overwritewithblock) {
                    md.write(path);// overwrite with active block only,
                } // other blocks were dismissed
                else {
                    md.writeBlock(path);// either if save active block or all, save active, other blocks where already managed
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }// function saveMd

    public MetaData getMetaData(String block) {
        return new MetaData(block + "@" + filename);
    }

    public boolean hasClasses() {
        return isClassification;
    }

    public boolean getNormalized() {
        return normalize;
    }

   

    public MetaData getMd(List<Long> ids) {
        MetaData selmd = null;
        try {
            long[] ids2 = new long[ids.size()];
            for (int i = 0; i < ids.size(); i++) {
                ids2[i] = ids.get(i);
            }
            selmd = new MetaData();
            selmd.importObjects(md, ids2);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return selmd;
    }

   
    public void setWindow(GalleryJFrame frame) {
        window = frame;
    }

    public void write(String path) {
        md.write(path);
    }

  

    public MetaData getMd() {
        return md;
    }

    public boolean isRenderLabel(ColumnInfo ci) {
        if(!renderImages)
            return false;
        if (renderLabel.equals("first") )
        {
            if(md instanceof ScipionMetaData)
                return ci.render && ci.visible;
            else
                return MetaData.isImage(ci.label) && ci.visible;
        }
        
        for (String i : getRenderLabels()) {
            if (i.equals(ci.labelName) && ci.visible) {
                return true;
            }
        }
        return false;
    }

    public boolean isVisibleLabel(ColumnInfo ci) {
        if (getVisibleLabels() == null) {
            return true;
        }
        for (String i : getVisibleLabels()) {
            if (i.equals(ci.labelName)) {
                return true;
            }
        }
        return false;
    }

    public void orderLabels() {
        String[] orderLabels = getOrderLabels();
        if (orderLabels == null) {
            return;
        }

        ColumnInfo aux;
        int j, k = 0;//k index added to avoid errors if some order columns are not present. This way order index is k not i
        for (int i = 0; i < orderLabels.length; i++) {
        	j = 0;
            for (ColumnInfo ci : labels) {
                if (ci.labelName.equals(orderLabels[i])) {

                    aux = labels.get(k);
                    labels.set(k, ci);
                    labels.set(j, aux);
                    k ++;
                }
                j ++;
            }
        }
    }

    public void removeCTF(int row) {
        md.removeCTF(ids[row]);
    }

    public boolean isRecalculateCTF(int row) {
        if(!isCTFMd())
            return false;
        long id = ids[row];
        if(md.containsCTF(id))
            return true;
        return false;
    }
    
    public boolean hasRecalculateCTF()
    {
        return md.hasRecalculateCTF();
    }
    
    public void recalculateCTF(int row, boolean[] selection, EllipseCTF ellipseCTF, String sortFn) 
    {
        if(isEnabled(row))
        {
            for(int i = 0; i <selection.length; i ++ )
                if(selection[i] && isEnabled(i))
                    md.putCTF(ids[i], ellipseCTF);
            EstimateFromCTFTask estimateFromCTFTask = new EstimateFromCTFTask(
                    ellipseCTF, 90 - ellipseCTF.getDefocusAngle(), 
                    md.getPSDFile(ids[row]), ellipseCTF.getD(), window.getTasksEngine(), row, sortFn);
            window.getTasksEngine().add(estimateFromCTFTask);
        }
    }
    
     public void showCTF(boolean profile, int row, boolean[] selection, TasksEngine ctfTasks) {
        try {
            long id = ids[row];
            
            String psd = md.getPSDFile(id);
            String psden = md.getPSDEnhanced(id);

            boolean noPsdEnhanced = (psden == null);
            // Use the same PSD image in case no PSD enhanced

            if (noPsdEnhanced)
                psden = psd;

            if (profile)
            {
                if (noPsdEnhanced)
                    new CTFAnalyzerJFrame(psden, md.getCTFDescription(id), psd, md.getEllipseCTF(id).getSamplingRate());
                else
                    new CTFAnalyzerJFrame(psden, md.getCTFFile(id), psd);
            }
            else
            {
                ImageGeneric img = new ImageGeneric(psden);
                ImagePlus imp = XmippImageConverter.readToImagePlus(img);
                img.destroy();
                EllipseCTF ctfparams = md.getEllipseCTF(id, imp.getWidth());
                String sortfn = createSortFile(psd, row);
                XmippUtil.showImageJ(Tool.VIEWER); // removed Toolbar.FREEROI
                CTFRecalculateImageWindow ctfiw = new CTFRecalculateImageWindow(this, selection, imp, psd, ctfparams,
                                                                                ctfTasks, row, sortfn);
            }

        } catch (Exception e) {
            XmippDialog.showError(window, e.getMessage());
        }
    }
     
    public Geometry getGeometry(long id, ColumnInfo ci)
    {
        return getGeometry(id, "2D", ci);
    }

    
    public Geometry getGeometry(long id, String type, ColumnInfo ci)
    {
        if(md.containsLabel(MDLabel.MDL_IMAGE) && ci.label != MDLabel.MDL_IMAGE)
        	return null;
        
        if(md.containsLabel(MDLabel.RLN_IMAGE_NAME) && ci.label != MDLabel.RLN_IMAGE_NAME)
        	return null;
        
        
        double shiftx, shifty, psiangle, scaleFactor=1;
        boolean flip;
        
        if(md.containsLabel(MetaData.GEOMETRY_LABELS))
        {
	        shiftx = md.getValueDouble(MDLabel.MDL_SHIFT_X, id);
	        shifty = md.getValueDouble(MDLabel.MDL_SHIFT_Y, id);
	        psiangle = md.getValueDouble(MDLabel.MDL_ANGLE_PSI, id);
	        flip = md.getValueBoolean(MDLabel.MDL_FLIP, id);
	        
	        if(md.containsLabel(MDLabel.MDL_SCALE))
	        	scaleFactor = md.getValueDouble(MDLabel.MDL_SCALE, id);
        }
        else if(md.containsLabel(MetaData.GEOMETRY_RELION_LABELS))
        {
	        shiftx = md.getValueDouble(MDLabel.RLN_ORIENT_ORIGIN_X, id);
	        shifty = md.getValueDouble(MDLabel.RLN_ORIENT_ORIGIN_Y, id);
	        psiangle = -1 * md.getValueDouble(MDLabel.RLN_ORIENT_PSI, id);
	        flip = md.getValueBoolean(MDLabel.MDL_FLIP, id);
        }
        else
        	return null;
        
        return new Geometry(shiftx, shifty, psiangle, flip, scaleFactor);
    }
        
	public void setRenderLabels(String[] renderLabels) {
            this.renderLabels = renderLabels;
    }

    public void setVisibleLabels(String[] visibleLabels) {
        this.visibleLabels = visibleLabels;
    }

    public void setOrderLabels(String[] orderLabels) {
        this.orderLabels = orderLabels;
    }

    public void setRenderLabel(String renderLabel) {
        this.renderLabel = renderLabel;
    }
    
    
    public String[] getRenderLabels()
    {
       return renderLabels;
    }
   
   public String[] getVisibleLabels()
   {
       return visibleLabels;
   }
   
   public String[] getOrderLabels()
   {
       return orderLabels;
   }
    
   public Long getSelVolId()
   {
       if(selectedVolFn == null || selectedVolFn.isEmpty())
           return null;
       String vol;
       for(Long id: ids)
       {
           vol = md.getValueString(ciFirstRender.label, id);
           if(vol.equals(selectedVolFn))
               return id;
       }
       return null;
   }
   
   public MetaData getImagesMd() {
	   return getImagesMd(null, false);
   }
   
   public MetaData getImagesMd(boolean[] selection, boolean selected) {
        int idlabel = getRenderLabel();
        MDRow mdRow;
        MetaData imagesmd = new MetaData();
        int index = 0;
        String imagepath;
        long imageid;
        for (long id : md.findObjects()) {
            if (isEnabled(index) && (!selected || selection[index])) {
                imagepath = md.getValueString(idlabel, id, true);
                if (imagepath != null && ImageGeneric.exists(imagepath)) {
                    imageid = imagesmd.addObject();
                    if (useGeo()) 
                    {
                    	Geometry geo = getGeometry(id, ciFirstRender);
                        mdRow = new MDRow();
                        //md.getRow(mdRow, id);//copy geo info in mdRow
                        mdRow.setValueDouble(MDLabel.MDL_SHIFT_X, geo.shiftx);
                        mdRow.setValueDouble(MDLabel.MDL_SHIFT_Y, geo.shifty);
                        mdRow.setValueDouble(MDLabel.MDL_ANGLE_PSI, geo.psiangle);
                        mdRow.setValueBoolean(MDLabel.MDL_FLIP, geo.flip);
                        imagesmd.setRow(mdRow, imageid);
                    }
                    imagesmd.setValueString(MDLabel.MDL_IMAGE, imagepath, imageid);
                }
            }
            index++;
        }
        return imagesmd;
   }
   
    
    
    public String[] getSortBy()
    {
        return sortby;
    }
    
    
    public boolean isChimeraClient()
    {
        return parameters.getChimeraPort() != null;
    }
    
    public boolean hasDisabled()
    {
    	for(int i = 0; i < ids.length; i ++)
    		if(!isEnabled(i))
    			return true;
    	return false;
    }

	public int getEnabledCount()
	{
		int count = 0;
		for(int i = 0; i < ids.length; i ++)
    		if(isEnabled(i))
    			count ++;
		return count;
	}
    
	public String getPreffix() {
        return selectedBlock;
    }
    
	public void setZoom(int zoom)
	{
		this.zoom = zoom;
	}
        
	public boolean allowsVolumeMode()
	{
		boolean allowsVolumeMode = !selectedVolFn.isEmpty();
		return allowsVolumeMode;
	}
}// class GalleryData
