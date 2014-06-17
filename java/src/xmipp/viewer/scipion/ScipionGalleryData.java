/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.Color;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.Params;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.models.ClassInfo;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryData extends GalleryData{

    public ScipionGalleryData(ScipionGalleryJFrame window, String fn, Params parameters, ScipionMetaData md) {
        super(window, fn, parameters, md);
        
        mdBlocks = md.getBlocks();
        selectedBlock = mdBlocks[0];
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

    
    public String getValueFromLabel(int index, int label)
    {
        return ((ScipionMetaData)md).getValueFromLabel(index, label);
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
    @Override
    public ScipionMetaData getSelectionMd() {
            ScipionMetaData selectionMd = null;
            try
            {
                selectionMd = ((ScipionMetaData)md).getMd(getSelObjects());
                return selectionMd;
            } catch (Exception e) {
                    e.printStackTrace();
            }
            return selectionMd;
    }
        
    	/** Get all the images assigned to all selected classes */
	public MetaData getClassesImages(){
		ScipionMetaData mdImages = null; 
		MetaData md;
		for (int i = 0; i < ids.length; ++i){
			if (selection[i] && isEnabled(i)){
				md = getClassImages(i);
                                
				if(md != null)
                                {
                                    if(mdImages == null)
                                        mdImages = ((ScipionMetaData)md).getStructure("");
                                    mdImages.unionAll(md);
                                    md.destroy();
                                }
			}
		}
                
		return mdImages;
	}   
        
         	/** Get all the images assigned to all selected classes */
	public MetaData getEnabledClassesImages(){
		ScipionMetaData mdImages = null; 
		MetaData md;
		for (int i = 0; i < ids.length; ++i){
			if (isEnabled(i)){
				md = getClassImages(i);
                                
				if(md != null)
                                {
                                    if(mdImages == null)
                                        mdImages = ((ScipionMetaData)md).getStructure("");
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
                        ScipionMetaData childmd = ((ScipionMetaData)md).getEMObject(id).childmd;
			return childmd.getMd(childmd.getEnabledObjects());
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
                                return String.format("Class %s (%d images)", emo.getId(), emo.childmd.size());
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
        
        
	/** Get the assigned class of some element */
	public ClassInfo getItemClassInfo(int index) {
		return null;
	}

	/** Set item class info in md */
	private void setItemClassInfo(long id, ClassInfo cli) {
		
	}

	/** Set the class of an element */
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

	/** Load classes structure if previously stored */
	public void loadClassesInfo() {
		
	}// function loadClassesInfo
        

        public MetaData[] getClassesMd() 
        {
            return null;
        }
        
        /** Add a new class */
	public void addClass(ClassInfo ci) {
		
	}

	/** Remove a class from the selection */
	public void removeClass(int classNumber) 
        {
		
	}
        
        public boolean hasClasses()//for Scipion usage only
        {
            return mdBlocks.length > 1 && ((ScipionMetaData)md).getSelf().contains("Class");
        }
         
        public boolean isCTFMd() {
		
		return false;//fixme?? provide ctf information maybe
	}
        
        public boolean hasMicrographParticles() {
		return false;//fixme?? cannot open picker from sqlite
	}
        
        public List<ScipionMetaData.EMObject> getSelObjects()
        {
            List<ScipionMetaData.EMObject> emos = new ArrayList<ScipionMetaData.EMObject>();
            for(int i = 0; i < selection.length; i ++)
                if(selection[i] && md.getEnabled(ids[i]))
                    emos.add(((ScipionMetaData)md).getEMObjects().get(i));
            return emos;
        }
        
        public List<ScipionMetaData.EMObject> getEnabledObjects()
        {
            return ((ScipionMetaData)md).getEnabledObjects();
        }
        
        /** This is only needed for metadata table galleries */
	public boolean isFile(ColumnInfo ci) {
		return ci.labelName.contains("filename");
	}
        
        public boolean isImageFile(ColumnInfo ci) {
            return ci.allowRender;
	}
        
        public MetaData getMd(List<Long> ids)
        {
            MetaData selmd = null;
            try {
                    long[] ids2 = new long[ids.size()];
                    for(int i = 0; i < ids.size(); i ++)
                        ids2[i] = ids.get(i);
                    selmd = ((ScipionMetaData)md).getStructure("");
                    selmd.importObjects(md, ids2);
                } catch (Exception e) {
                        e.printStackTrace();
                }
            return selmd;
        }
        
        @Override
        public void openMetadata(MetaData md)
        {
            ScipionGalleryData data2 = new ScipionGalleryData(null, md.getFilename(), new ScipionParams(), (ScipionMetaData)md);
            ScipionGalleryJFrame frame = new ScipionGalleryJFrame(data2);
            data2.setWindow(frame);
        }
        
        
       
        public void overwrite(String path)
        {
            ((ScipionMetaData)md).overwrite(filename, path);
        }

    public String getScipionType() {
        if(hasClasses())
            return "Particle";
        String self = ((ScipionMetaData)md).getSelf();
        if(self.equals("CTFModel"))
            return "Micrograph";
        return self;
    }
    
    public String getSelf() {
        
        return ((ScipionMetaData)md).getSelf();
    }

    public String getPrefix() {
        return ((ScipionMetaData)md).getPrefix();
    }
    
}
