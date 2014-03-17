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
import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippStringUtils;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryJFrame extends GalleryJFrame {

    private final String type;
    private final String script;
    private final String projectid;
    private JButton cmdbutton;
    private String selectionmdfile;
    private JButton classcmdbutton;
    private String firebrick = "#B22222";
    private String lightgrey = "#EAEBEC";
    private String python;
    private String inputimagesid;
    private String inputid;

    public ScipionGalleryJFrame(String filename, MetaData md, ScipionParams parameters) {
        super(filename, md, parameters);
        try {
            type = parameters.type;
            python = parameters.python;
            script = parameters.script;
            projectid = parameters.projectid;
            inputid = parameters.inputid;
            inputimagesid = parameters.inputimagesid;
            selectionmdfile = String.format("%s%sselection%s", projectid, File.separator, getFileExtension());
            
            initComponents();
        } catch (Exception ex) {
            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
            throw new IllegalArgumentException(ex.getMessage());
        }
    }

    private void initComponents() {
        if (type != null) {
            final String inputType = is2DClassificationMd()? "Classes2D": type;
            cmdbutton = getScipionButton("Create New Set Of " + type, new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    try {
                        
                        if(is2DClassSelection())
                            saveImagesFromClassSelection(selectionmdfile);
                        else
                            saveSelection(selectionmdfile, true);
                        String[] command = new String[]{python, script, selectionmdfile, inputType, type, projectid, inputid, inputimagesid};
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
                            String[] command = new String[]{python, script, selectionmdfile, "Classes2D", "Classes2D", projectid, inputid, inputimagesid};
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
    
    protected void runCommand(final String[] command)
    {
        XmippWindowUtil.blockGUI(ScipionGalleryJFrame.this, "Creating set ...");
        new Thread(new Runnable() {

            @Override
            public void run() {

                try {
                    
                    String output = XmippUtil.executeCommand(command);
                    XmippWindowUtil.releaseGUI(ScipionGalleryJFrame.this.getRootPane());
                    if(output != null && !output.isEmpty())
                    {
                        XmippDialog.showInfo(ScipionGalleryJFrame.this, output);
                        System.out.println(output);
                    }
                    
                } catch (Exception ex) {
                    throw new IllegalArgumentException(ex.getMessage());
                }

            }
        }).start();
    }
    
    public JButton getScipionButton(String text, ActionListener listener)
    {
    
        JButton button = new JButton(text);
        button.addActionListener(listener);
   
       
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
        Color forecolor = isenabled? Color.WHITE: Color.GRAY;
        cmdbutton.setEnabled(isenabled);
        cmdbutton.setBackground(color);
        cmdbutton.setForeground(forecolor);
        if(classcmdbutton != null)
        {
            isenabled = is2DClassificationMd() && isenabled;
            color = Color.decode(isenabled? firebrick: lightgrey); 
            forecolor = isenabled? Color.WHITE: Color.GRAY;
            classcmdbutton.setEnabled( isenabled);
            classcmdbutton.setBackground(color);
            classcmdbutton.setForeground(forecolor);
        }
    }

}
