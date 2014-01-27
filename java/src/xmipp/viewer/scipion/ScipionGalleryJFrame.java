/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import xmipp.jni.MetaData;
import xmipp.utils.Param;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryJFrame extends GalleryJFrame {
    private final String cmdname;
    private final String cmdscript;

    public ScipionGalleryJFrame(String filename, MetaData md, Param parameters) {
        super(filename, md, parameters);
        cmdname = parameters.cmdname;
        cmdscript = parameters.cmdscript;
        initComponents();
    }

    private void initComponents() {
        if(cmdname != null)
            commandspn.add(XmippWindowUtil.getTextButton(cmdname, new ActionListener() {

                       @Override
                       public void actionPerformed(ActionEvent ae) {
                           try {
                               Runtime.getRuntime().exec(cmdscript);
                           } catch (IOException ex) {
                               
                               Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
                           }
                       }
                   }));
    }
    
   
}
