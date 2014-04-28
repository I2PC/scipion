/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.Window;
import java.util.ArrayList;
import xmipp.jni.Filename;
import xmipp.jni.MetaData;
import xmipp.utils.Params;
import xmipp.viewer.models.ColumnInfo;
import xmipp.viewer.models.GalleryData;

/**
 *
 * @author airen
 */
public class ScipionGalleryData extends GalleryData{

    public ScipionGalleryData(Window window, String fn, Params parameters, MetaData md) {
        super(window, fn, parameters, md);
        selectedBlock = ((ScipionMetaData)md).getSelf() + "s";
        mdBlocks = new String[]{selectedBlock};
    }
    
    public void setFileName(String file) {
		filename = file;
                


    }
    
    public void loadLabels() {
        
            ScipionMetaData smd = (ScipionMetaData)md;
            if(zoom == 0)
                zoom = 100;
            labels = smd.getColumnsInfo();
            orderLabels();
            for(ColumnInfo ci :labels)
                if(ci.render && ciFirstRender == null)
                {
                    ciFirstRender = ci;
                    break;
                }
            
    
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
    
}
