/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.MetaData;
import xmipp.utils.Param;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryJFrame extends GalleryJFrame {

    private final String type;
    private final String script;
    private final String projectid;
    private final String imagesid;
    private JButton cmdbutton;
    private String selectionmdfile;

    public ScipionGalleryJFrame(String filename, MetaData md, ScipionParams parameters) {
        super(filename, md, parameters);
        type = parameters.type;
        script = parameters.script;
        projectid = parameters.projectid;
        imagesid = parameters.imagesid;
        selectionmdfile = "selection.xmd";
        selectionmdfile = new File(selectionmdfile).getAbsolutePath();
        initComponents();
    }

    private void initComponents() {
        if (type != null) {
            cmdbutton = XmippWindowUtil.getTextButton("Create New Set Of " + type, new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    try {
                        
                        if(is2DClassSelection())
                            saveImagesFromClassSelection(selectionmdfile);
                        else
                            saveSelection(selectionmdfile, true);
                        
                        XmippWindowUtil.blockGUI(ScipionGalleryJFrame.this, "Creating set ..." + "(You may need to refresh the main window to visualize output)");
                        new Thread(new Runnable() {

                            @Override
                            public void run() {
                                String command = String.format("%s %s %s %s %s", script, selectionmdfile, type, projectid, imagesid);
                                try {
                                    XmippUtil.executeCommand(command);
                                    XmippWindowUtil.releaseGUI(ScipionGalleryJFrame.this.getRootPane());
                                } catch (Exception ex) {
                                    throw new IllegalArgumentException(ex.getMessage());
                                }
                                
                            }
                        }).start();
                        
                    } catch (Exception ex) {
                        Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
                        
                    }
                }
            });

            buttonspn.add(cmdbutton);
            cmdbutton.setEnabled(false);
        }
       
    }
    
   

    public void selectItem(int row, int col)
    {
        super.selectItem(row, col);
        cmdbutton.setEnabled(isImageSelected());
    }

    protected void tableMouseClicked(MouseEvent evt)
    {
        super.tableMouseClicked(evt);
        cmdbutton.setEnabled(isImageSelected());
    }


}
