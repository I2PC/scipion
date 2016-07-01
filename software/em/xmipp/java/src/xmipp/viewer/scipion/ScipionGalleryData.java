/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import xmipp.utils.ScipionParams;

import java.io.File;
import java.sql.SQLException;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import xmipp.ij.commons.Geometry;
import xmipp.jni.EllipseCTF;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.Params;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.models.ClassInfo;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;
import xmipp.viewer.scipion.ScipionMetaData.EMObject;

/**
 *
 * @author airen
 */
public class ScipionGalleryData extends GalleryData {
    
    public ScipionGalleryData(ScipionGalleryJFrame window, String fn, Params parameters) {
        this(window, parameters, new ScipionMetaData(fn));
    }

    public ScipionGalleryData(ScipionGalleryJFrame window, Params parameters, ScipionMetaData md) {
        super(window, parameters, md);
    }

    public void setFileName(String file) {
        if (file.contains("@")) {
            int sep = file.lastIndexOf("@");
            selectedBlock = file.substring(0, sep);
            filename = file.substring(sep + 1);
        }
        filename = file;
        mdBlocks = ((ScipionMetaData)md).getBlocks();
        selectedBlock = mdBlocks[0];

    }

    
    @Override
    public ColumnInfo initColumnInfo(int label)
    {
        return ((ScipionMetaData)md).getColumnInfo(label);
    }

    public String getValueFromLabel(int index, int label) {
        return ((ScipionMetaData) md).getValueFromLabel(index, label);
    }

    
    

    /**
     * Create a metadata just with selected items
     */
    @Override
    public ScipionMetaData getSelectionMd(boolean[] selection) {
        return null;//metadata operations are not used in scipion
    }

    /**
     * Get all the images assigned to all selected classes
     */
    public MetaData getClassesImages() {
        return null;
    }

    

    /**
     * Get the metadata with assigned images to this classes
     */
    public MetaData getClassImages(int index) {
        try {
            long id = ids[index];
            ScipionMetaData childmd = ((ScipionMetaData) md).getEMObject(id).childmd;
            return childmd;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    public String getLabel(long objId, int label) {
        try {
            if (isClassification) {
                ScipionMetaData.EMObject emo = ((ScipionMetaData) md).getEMObject(objId);
                return String.format("Class %s (%d images)", emo.getId(), emo.childmd.size());
            } else {
                return md.getValueString(label, objId);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    public void readMd() {
        hasMdChanges = false;
        hasClassesChanges = false;
        md = getMetaData(selectedBlock);
    }

    public MetaData getMetaData(String block) {
        if (md.getBlock().equals(block)) {
            return md;
        }
        ScipionMetaData child = ((ScipionMetaData) md).getChild(block);
        if (child != null) {
            return child;
        }
        ScipionMetaData parent = ((ScipionMetaData) md).getParent();
        if (parent.getBlock().equals(selectedBlock))// from child to parent
        {
            return parent;
        }
        return parent.getChild(selectedBlock);

    }

    public ScipionMetaData getScipionMetaData(){
        return (ScipionMetaData) md;
    }

    /**
     * Get the assigned class of some element
     */
    public ClassInfo getItemClassInfo(int index) {
        return null;
    }

    /**
     * Set item class info in md
     */
    private void setItemClassInfo(long id, ClassInfo cli) {

    }

    /**
     * Set the class of an element
     */
    public void setItemClass(int index, ClassInfo cli) {

    }

    public ClassInfo getClassInfo(int classNumber) {
        return null;
    }

    /**
     * Compute and update the number of classes and images assigned to this
     * superclass
     */
    public void updateClassesInfo() {

    }// function upateClassesInfo

    /**
     * Load classes structure if previously stored
     */
    public void loadClassesInfo() {

    }// function loadClassesInfo

    public MetaData[] getClassesMd() {
        return null;
    }

    /**
     * Add a new class
     */
    public void addClass(ClassInfo ci) {

    }

    /**
     * Remove a class from the selection
     */
    public void removeClass(int classNumber) {

    }

    public boolean hasClasses()//for Scipion usage only
    {
        return mdBlocks.length > 1 && ((ScipionMetaData) md).getSetType().contains("SetOfClasses");
    }

    public boolean hasMicrographParticles() {
        return false;//fixme?? cannot open picker from sqlite
    }

    public List<ScipionMetaData.EMObject> getEMObjects() {
        return ((ScipionMetaData)md).getEMObjects();
    }

    /**
     * This is only needed for metadata table galleries
     */
    public boolean isFile(ColumnInfo ci) {
        return ci.labelName.contains("filename") || ci.labelName.equals("_psdFile");
    }

    public boolean isImageFile(ColumnInfo ci) {
        return ci.allowRender;
    }


    public void overwrite(String path, boolean[] selection) throws SQLException {
        ((ScipionMetaData) md).overwrite(filename, path, selection);
    }

    

    public String getPreffix() {
        return ((ScipionMetaData) md).getPreffix();
    }
    
    @Override
    public void removeCTF(int row) {
        ScipionMetaData.EMObject emo = ((ScipionMetaData) md).getEMObjects().get(row);
        emo.setLabel("");
        emo.setComment("");
        super.removeCTF(row);
        window.fireTableRowsUpdated(row, row);
    }

    public String createSortFile(String psdFile, int row) {
        return null;
    }
    
    @Override
    public void recalculateCTF(int row, boolean[] selection, EllipseCTF ellipseCTF, String sortFn) {
        
        ScipionMetaData.EMObject emo;
        long id;
        if(isEnabled(row))
        {
            for(int i = 0; i <selection.length; i ++ )
                if(selection[i] && isEnabled(i))
                {
                    id = ids[i];
                    emo = ((ScipionMetaData) md).getEMObject(id);
                    md.putCTF(id, ellipseCTF);
                    emo.setLabel("(recalculate ctf)");
                    emo.setComment(ellipseCTF.toString());
                    window.fireTableRowsUpdated(i, i);
                }
        }
        
    }
    
    public ColumnInfo getGeoMatrixColumn(ColumnInfo ci)
    {
        return ((ScipionMetaData)md).getGeoMatrixColumn(ci);
    }
    
    public Geometry getGeometry(long id, ColumnInfo ci)
    {
        return getGeometry(id, "2D", ci);
    }
            
    public Geometry getGeometry(long id, String type, ColumnInfo ci)
    {
        if (!containsGeometryInfo(type)) 
            return null;
        ScipionMetaData.EMObject emo = ((ScipionMetaData)md).getEMObject(id);
        ColumnInfo column = getGeoMatrixColumn(ci);
        
        String matrix = (String)emo.getValue(column);
        return new Geometry(matrix);
    }

    public int getEnabledCount() {
        return ((ScipionMetaData)md).getEnabledCount();
    }
        
     /**
     * Set enabled state
     */
    @Override
    public void setEnabled(int index, boolean isenabled) {
        try {
            if (!isVolumeMode()) { // slices in a volume are always enabled
                getEMObjects().get(index).setEnabled(isenabled);
                hasMdChanges = true;
                if(!isenabled && isRecalculateCTF(index))
                    removeCTF(index);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public boolean hasMdChanges()
    {
        return ((ScipionMetaData)md).isChanged();
    }
    
    /**
     * Check if an item is enabled or not
     */
    public boolean isEnabled(int index) {
        try {
            if (isVolumeMode()) {
                return true;
            }
            return getEMObjects().get(index).isEnabled();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return true;
    }
    
    public String getValueFromCol(int index, ColumnInfo ci) {
        try {
            return getEMObjects().get(index).getValueString(ci);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

     public void setValueToCol(int index, ColumnInfo ci, String value) {
        try {
            getEMObjects().get(index).setValue(ci, value);
            setMdChanges(true);
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void loadSelection(String selectedPath) {
        ((ScipionMetaData)md).loadSelection(selectedPath);
            
    }
    
       
    public String getExtraPath()
    {
        if(filename == null)
            return null;
        String dir = new File(filename).getParent();
        if(dir == null)
            return null;
        return Filename.join(dir, "extra");
    }
    
        /**
     * Return true if current metadata comes from 2d classification
     */
    public boolean checkifIsClassificationMd() {
        return ((ScipionMetaData)md).isClassificationMd();
    }
    
    public void runObjectCommand(int index, String objectCommand) {
        try {
            ScipionParams params = (ScipionParams)parameters;
            String command = String.format("run function runObjectCommand '%s' %s %s", objectCommand, params.inputid, getId(index));
            XmippWindowUtil.runCommand(command, params.port);
//            String[] cmd = new String[]{params.python, params.getObjectCmdScript(), String.format("'%s'", objectCommand), params.projectid, params.other, String.valueOf(getId(index))};
//            XmippWindowUtil.executeCommand(cmd, false);
        } catch (Exception ex) {
            Logger.getLogger(GalleryData.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
    
    
    
    public MetaData getImagesMd(boolean[] selection, boolean selected) {
    		
            MetaData imagesmd = new MetaData();
            int index = 0;
            String imagepath;
            EMObject emo;
            long imageid;
            String matrix;
            
            for (long id : md.findObjects()) {
                if (isEnabled(index) && (!selected ||selection[index])) {
                        
                    imagepath = md.getValueString(getRenderLabel(), id, true);
                    if (imagepath != null && ImageGeneric.exists(imagepath)) {
                        imageid = imagesmd.addObject();
                        imagesmd.setValueString(MDLabel.MDL_IMAGE, imagepath, imageid);
                        if (useGeo()) 
                        {
                            emo = ((ScipionMetaData)md).getEMObject(id);
                            matrix = String.format("'%s'", emo.getValueString(getGeoMatrixColumn(ciFirstRender)));
                            imagesmd.setValueString(MDLabel.MDL_TRANSFORM_MATRIX, matrix, imageid);//copy geo info in mdRow
                        }
                        
                    }
                }
                index++;
            }
            return imagesmd;
    }
    
    public String getChimeraProjectionCmd(int row) {
    	int label = getRenderLabel();
		ColumnInfo ci = getColumnInfo(label);
        Geometry geo = getGeometry(ids[row], "3D", ci);
        if(geo == null)
            geo = getGeometry(ids[row], "Projection", ci);
        String command = String.format("rotate_matrix '%s'", geo.getMatrix());
        return command;
    }

	
    
    
    
}
