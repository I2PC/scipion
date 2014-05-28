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
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryJFrame extends GalleryJFrame {

    private String type;
    private String script;
    private String projectid;
    private JButton cmdbutton;
    private String selectionmdfile;
    private JButton classcmdbutton;

    private String python;
    private String inputimagesid;
    private String inputid;
    private HashMap<String, String> msgfields;
    private final String runNameKey = "Run name:";
    
   

    public ScipionGalleryJFrame(String filename, MetaData md, ScipionParams parameters) {
        super(filename, md, parameters);
        readScipionParams(parameters);
       
    }
    
      public ScipionGalleryJFrame(ScipionGalleryData data) {
        super(data);
        readScipionParams((ScipionParams)data.parameters);
    }

    protected void readScipionParams(ScipionParams parameters)
    {
        try {
            type = parameters.type;
            python = parameters.python;
            script = parameters.script;
            projectid = parameters.projectid;
            inputid = parameters.inputid;
            inputimagesid = parameters.inputimagesid;
            selectionmdfile = String.format("%s%sselection%s", projectid, File.separator, getFileExtension());
            msgfields = new HashMap<String, String>();
            msgfields.put(runNameKey, "ProtUserSubset");
            

            initComponents();
        } catch (Exception ex) {
            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
            throw new IllegalArgumentException(ex.getMessage());
        }
    }
    private void initComponents() {
        if (type != null) {
            boolean isclassificationmd = isClassificationMd();
            String output = isclassificationmd? "Particles": type;   
            cmdbutton = getScipionButton("Create " + output, new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    saveSelection();
                    MetaData selectionmd = new MetaData(selectionmdfile);
                    String question = String.format("<html>Are you sure you want to create a new set of %s with <font color=red>%s</font> %s?", type, selectionmd.size(), (selectionmd.size() > 1)?"elements":"element");
                    ScipionMessageDialog dlg = new ScipionMessageDialog(ScipionGalleryJFrame.this, "Question", question, msgfields);
                    int create = dlg.action;
                    if (create == ScipionMessageDialog.OK_OPTION) {
                        createSubset(dlg.getFieldValue(runNameKey));
                    }
                    selectionmd.destroy();
                }
            });

            if(isclassificationmd)
            {
                classcmdbutton = getScipionButton("Create Classes", new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        saveClassSelection(selectionmdfile);
                        MetaData selectionmd = new MetaData(selectionmdfile);
                        String msg = String.format("<html>Are you sure you want to create a new set of Classes with <font color=red>%s</font> %s?", selectionmd.size(), (selectionmd.size() > 1)?"elements":"element");
                        ScipionMessageDialog dlg = new ScipionMessageDialog(ScipionGalleryJFrame.this, "Question", msg, msgfields);
                        int create = dlg.action;
                        if (create == ScipionMessageDialog.OK_OPTION) {
                            createSubsetOfClasses(dlg.getFieldValue(runNameKey));
                        }
                        selectionmd.destroy();

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

    protected void runCommand(final String[] command) {
        XmippWindowUtil.blockGUI(ScipionGalleryJFrame.this, "Creating set ...");
        new Thread(new Runnable() {

            @Override
            public void run() {

                try {

                    String output = XmippUtil.executeCommand(command);
                    XmippWindowUtil.releaseGUI(ScipionGalleryJFrame.this.getRootPane());
                    if (output != null && !output.isEmpty()) {
                        XmippDialog.showInfo(ScipionGalleryJFrame.this, output);
                        System.out.println(output);
                    }

                } catch (Exception ex) {
                    throw new IllegalArgumentException(ex.getMessage());
                }

            }
        }).start();
    }

    public JButton getScipionButton(String text, ActionListener listener) {

        JButton button = new JButton(text);
        button.addActionListener(listener);

        return button;
    }

    public void selectItem(int row, int col) {
        super.selectItem(row, col);
        enableActions();

    }

    protected void tableMouseClicked(MouseEvent evt) {
        super.tableMouseClicked(evt);
        enableActions();
    }

    protected void enableActions() {
        boolean isenabled = isImageSelected();
        Color color = Color.decode(isenabled ? ScipionMessageDialog.firebrick : ScipionMessageDialog.lightgrey);
        Color forecolor = isenabled ? Color.WHITE : Color.GRAY;
        if(cmdbutton != null)
        {
            cmdbutton.setEnabled(isenabled);
            cmdbutton.setBackground(color);
            cmdbutton.setForeground(forecolor);
        }
        if(classcmdbutton != null)
        {
            isenabled = isClassificationMd() && isenabled;
            color = Color.decode(isenabled? ScipionMessageDialog.firebrick: ScipionMessageDialog.lightgrey); 
            forecolor = isenabled? Color.WHITE: Color.GRAY;
            classcmdbutton.setEnabled( isenabled);

            classcmdbutton.setBackground(color);
            classcmdbutton.setForeground(forecolor);
        }
    }

    public void createSubset(String runname) {
        try {
            String output = type;     
            if(isClassSelection())
            {
                saveImagesFromClassSelection(selectionmdfile);
                output = "Particles";
            }
            else
                saveSelection(selectionmdfile, true);
            
            String[] command = new String[]{python, script, runname, selectionmdfile, type, output, projectid, inputid, inputimagesid};

            runCommand(command);

        } catch (Exception ex) {
            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);

        }
    }
    
    public void saveSelection() {
        try {
            if (isClassSelection()) {
                saveImagesFromClassSelection(selectionmdfile);
            } else {
                saveSelection(selectionmdfile, true);
            }
        } catch (Exception ex) {
            throw new IllegalArgumentException(ex.getMessage());
        }
    }

    
    public void createSubsetOfClasses(String runname)
    {
        try 
        {

            saveClassSelection(selectionmdfile);
            String[] command = new String[]{python, script, runname, selectionmdfile, type, type, projectid, inputid, inputimagesid};

            runCommand(command);
        } catch (Exception ex) {
            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    

}
