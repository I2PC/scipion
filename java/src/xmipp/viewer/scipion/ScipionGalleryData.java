/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.util.logging.Level;
import java.util.logging.Logger;
import xmipp.jni.MetaData;
import xmipp.utils.Params;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;

/**
 *
 * @author airen
 */
public class ScipionGalleryData extends GalleryData{

    public ScipionGalleryData(ScipionGalleryJFrame window, String fn, Params parameters, ScipionMetaData md) {
        super(window, fn, parameters, md);
        
        mdBlocks = md.getBlocks();
        selectedBlock = mdBlocks[0];
        isClassification = mdBlocks.length > 1;//FIXME:Temporarily
        classes = ((ScipionMetaData)md).getClasses();
    }
    
    public void setFileName(String file) {
        if(file.contains("@"))
        {
            int sep = file.lastIndexOf("@");
            selectedBlock = file.substring(0, sep);
            filename = file.substring(sep + 1);
        }
	filename = file;
                
    }
    
    public void loadLabels() {
            ScipionMetaData smd = (ScipionMetaData)md;
            
            labels = smd.getColumnsInfo();
            ciFirstRender = null;
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
        getSelectionMd().write(path);
        
    }
    
     public void saveClassSelection(String path)
     {
        try {
            saveSelection(path, true);//Scipion metadata saves recursively
        } catch (Exception ex) {
            Logger.getLogger(ScipionGalleryData.class.getName()).log(Level.SEVERE, null, ex);
        }
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
                selectionMd = ((ScipionMetaData)md).getSelectionMd(getSelIds());
                return selectionMd;
            } catch (Exception e) {
                    e.printStackTrace();
            }
            return selectionMd;
    }
        
    	/** Get all the images assigned to all selected classes */
	public MetaData getSelClassesImages(){
		ScipionMetaData mdImages = null; 
		MetaData md;
		for (int i = 0; i < ids.length; ++i){
			if (selection[i] && isEnabled(i)){
				md = getClassImages(i);
                                
				if(md != null)
                                {
                                    if(mdImages == null)
                                        mdImages = ((ScipionMetaData)md).getStructure("Classes", "Objects");
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
        
        public String getLabel(long objId, int label)
        {
            try
                    {
                            if (isClassification)
                            {
                                ScipionMetaData.EMObject emo = ((ScipionMetaData)md).getEMObject(objId);
                                return String.format("Class %s (%d images)", emo.id, emo.childmd.size());
                            }
                            else
                                    return md.getValueString(label, objId);
                    }
                    catch (Exception e)
                    {
                            e.printStackTrace();
                    }
                    return null;
        }
        
        public void readMd() {
                hasMdChanges = false;
                hasClassesChanges = false;
                md = getMetaData(selectedBlock);
	}

        public MetaData getMetaData(String block)
        {
            if(md.getBlock().equals(block))
                    return md;
                ScipionMetaData child = ((ScipionMetaData)md).getChild(block);
                if(child != null)
                    return child;
                ScipionMetaData parent = ((ScipionMetaData)md).getParent();
                if(parent.getBlock().equals(selectedBlock))// from child to parent
                    return parent;
                return parent.getChild(selectedBlock);
                
        }

}
