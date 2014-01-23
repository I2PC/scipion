/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package scipion.viewer;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import xmipp.jni.MetaData;
import xmipp.utils.Param;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryJFrame extends GalleryJFrame {

    public ScipionGalleryJFrame(String filename, MetaData md, Param parameters) {
        super(filename, md, parameters);
        initComponents();
    }

    private void initComponents() {
         commandspn.add(XmippWindowUtil.getTextButton("Command", new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
                    }
                }));
    }
    
   
}
