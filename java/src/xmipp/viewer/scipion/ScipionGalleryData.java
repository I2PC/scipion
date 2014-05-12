/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.Color;
import java.util.List;
import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.Params;
import xmipp.viewer.models.ClassInfo;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;

/**
 *
 * @author airen
 */
public class ScipionGalleryData extends GalleryData{

    public ScipionGalleryData(ScipionGalleryJFrame window, String fn, Params parameters, ScipionMetaData md) {
        super(window, fn, parameters, md);
        
        classes = ((ScipionMetaData)md).getClasses();
        selectedBlock = classes[0].getComment();
        mdBlocks = new String[]{selectedBlock};
        isClassification = mdBlocks.length > 1;//FIXME:Temporarily
    }
    
    public void setFileName(String file) {
	filename = file;
                
    }
    
    public void loadLabels() {
        
            ScipionMetaData smd = (ScipionMetaData)md;
            
            labels = smd.getColumnsInfo();
            orderLabels();
            for(ColumnInfo ci :labels)
            {
                
                ci.visible = isVisibleLabel(ci);
                ci.render = isRenderLabel(ci);
                if(ci.render && ciFirstRender == null)
                    ciFirstRender = ci;
            }
            if(ciFirstRender == null)
                zoom = 100;
    }

    void setWindow(ScipionGalleryJFrame frame) {
        window = frame;
    }
    
    public String getValueFromLabel(int index, int label)
    {
        return ((ScipionMetaData)md).getValueFromLabel(index, label);
        
        
    }
    
    public boolean isFile(ColumnInfo ci) {
        return false;
    }
    
            /** Save selected items as a metadata */
    public void saveSelection(String path, boolean overwrite) throws Exception
    {
            ((ScipionMetaData)md).writeBlock(selectedBlock, path, getSelIds());
        
    }

    public boolean isColumnFormat()
    {
        return true;
    }
    

    /** Create a metadata just with selected items */
    public ScipionMetaData getSelectionMd() {
            ScipionMetaData selectionMd = null;
            try
            {
                selectionMd = new ScipionMetaData("selection.sqlite");
                selectionMd.importObjects(md, getSelIds());//import objects, if classes afterwords should import particles
            } catch (Exception e) {
                    e.printStackTrace();
            }
            return selectionMd;
    }
        
    	/** Get all the images assigned to all selected classes */
	public MetaData getSelClassesImages(){
		ScipionMetaData mdImages = new ScipionMetaData();
		MetaData md;
		for (int i = 0; i < ids.length; ++i){
			if (selection[i]){
				md = getClassImages(i);
				if(md != null)
                                {
                                    mdImages.unionAll(md);
                                    md.destroy();
                                }
			}
		}
		return mdImages;
	}   
        
        /** Get the metadata with assigned images to this classes */
	public MetaData getClassImages(int index) {
		try {
			long id = ids[index];
			return ((ScipionMetaData)md).getEMObject(id).childmd;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
        
       

        

}
