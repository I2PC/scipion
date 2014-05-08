/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.Window;
import java.io.File;
import xmipp.jni.MetaData;
import xmipp.utils.Params;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;
import xmipp.viewer.windows.SaveJDialog;

/**
 *
 * @author airen
 */
public class ScipionGalleryData extends GalleryData{

    public ScipionGalleryData(ScipionGalleryJFrame window, String fn, Params parameters, MetaData md) {
        super(window, fn, parameters, md);
        selectedBlock = ((ScipionMetaData)md).getSelf() + "s";
        mdBlocks = new String[]{selectedBlock};
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
	public void saveSelection() throws Exception
	{
		MetaData md = getSelectionMd();
		SaveJDialog dlg = new SaveJDialog(window, "selection.xmd", true);
		boolean save = dlg.showDialog();
		if (save)
		{
			boolean overwrite= dlg.isOverwrite();
			String path = dlg.getMdFilename();
			saveSelection(path, overwrite);
		}
		md.destroy();
	}
        
        	/** Save selected items as a metadata */
	public void saveSelection(String path, boolean overwrite) throws Exception
	{
		MetaData md = getSelectionMd();

                String file = path.substring(path.lastIndexOf("@") + 1, path.length());
                if (!new File(file).exists())// overwrite or append, save selection
                        md.write(path);
                else
                {
                        if (overwrite)
                                md.write(path);// overwrite with active block only, other
                                                                // blocks were dismissed
                        else
                                md.writeBlock(path);// append selection

                }
		
		md.destroy();
	}
        
        public boolean isColumnFormat()
        {
            return true;
        }
    
        
}
