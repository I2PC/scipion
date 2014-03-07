/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.UIManager;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.MetaData;
import xmipp.utils.Param;
import xmipp.utils.XmippDialog;
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
    private JButton classcmdbutton;
    private String firebrick = "#B22222";
    private String lightgrey = "#EAEBEC";

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
            cmdbutton = getScipionButton("Create New Set Of " + type, new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    try {
                        
                        if(is2DClassSelection())
                            saveImagesFromClassSelection(selectionmdfile);
                        else
                            saveSelection(selectionmdfile, true);
                        String command = String.format("%s %s %s %s %s", script, selectionmdfile, type, projectid, imagesid);
                        runCommand(command);
                        
                        
                    } catch (Exception ex) {
                        Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
                        
                    }
                }
            });
            if(is2DClassificationMd())
            {
                classcmdbutton = getScipionButton("Create New Set Of Classes", new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        try {
                            saveClassSelection(selectionmdfile);
                            String command = String.format("%s %s %s %s %s", script, selectionmdfile, "Classes2D", projectid, imagesid);
                            runCommand(command);
                        } catch (Exception ex) {
                            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                });
                buttonspn.add(classcmdbutton);
                
            }
            buttonspn.add(cmdbutton);
            enableActions();
            jcbBlocks.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    enableActions();
                }
            });
        }
       
    }
    
    protected void runCommand(final String command)
    {
        XmippWindowUtil.blockGUI(ScipionGalleryJFrame.this, "Creating set ...");
        new Thread(new Runnable() {

            @Override
            public void run() {

                try {
                    String output = XmippUtil.executeCommand(command);
                    XmippWindowUtil.releaseGUI(ScipionGalleryJFrame.this.getRootPane());
                    
                } catch (Exception ex) {
                    throw new IllegalArgumentException(ex.getMessage());
                }

            }
        }).start();
    }
    
    public JButton getScipionButton(String text, ActionListener listener)
    {
         
        //UIManager.getDefaults().put("Button.background", Color.decode(firebrick));
        UIManager.getDefaults().put("Button.foreground",Color.WHITE);
        //UIManager.getDefaults().put("Button.disabledShadow", Color.WHITE);
        
        JButton button = new JButton(text);
        
        button.addActionListener(listener);
        //button.setBackground(Color.decode(firebrick));
        button.setForeground(Color.WHITE);
       
        return button;
    }
    
   

    public void selectItem(int row, int col)
    {
        super.selectItem(row, col);
        enableActions();
       
    }

    protected void tableMouseClicked(MouseEvent evt)
    {
        super.tableMouseClicked(evt);
        enableActions();
    }

    protected void enableActions()
    {
        boolean isenabled = isImageSelected();
        Color color = Color.decode(isenabled? firebrick: lightgrey); 
        cmdbutton.setEnabled(isenabled);
        cmdbutton.setBackground(color);
        isenabled = is2DClassificationMd() && isenabled;
        color = Color.decode(isenabled? firebrick: lightgrey); 
        classcmdbutton.setEnabled( isenabled);
        classcmdbutton.setBackground(color);
         
    }

}
