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
    private String sqlitefile;
    private JButton classcmdbutton;
    private String python;
    private String inputid;
    private HashMap<String, String> msgfields;
    private final String runNameKey = "Run name:";
    private String other;
    
   

    public ScipionGalleryJFrame(String filename, ScipionMetaData md, ScipionParams parameters) {
        super(filename, md, parameters);
        readScipionParams(parameters);
        data.setWindow(this);
    }
    
      public ScipionGalleryJFrame(ScipionGalleryData data) {
        super(data);
        readScipionParams((ScipionParams)data.parameters);
        data.setWindow(this);
    }

    protected void readScipionParams(ScipionParams parameters)
    {
        try {
            type = ((ScipionGalleryData)data).getScipionType() + "s";
            python = parameters.python;
            script = parameters.script;
            projectid = parameters.projectid;
            inputid = parameters.inputid;
            sqlitefile = String.format("%s%sselection%s", projectid, File.separator, data.getFileExtension());
            msgfields = new HashMap<String, String>();
            msgfields.put(runNameKey, "ProtUserSubset");
            other = parameters.other;
            initComponents();
        } catch (Exception ex) {
            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
            throw new IllegalArgumentException(ex.getMessage());
        }
    }
    private void initComponents() {
        JButton closebt = getScipionButton("Close", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                close();
            }
        });
        buttonspn.add(closebt);
        if (type != null) {
            cmdbutton = getScipionButton("Create " + type, new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    int size;
                    MetaData imagesmd = null;    
                    if(data.hasClasses())
                        imagesmd = data.getEnabledClassesImages();
                    else
                        imagesmd = data.getMd(data.getEnabledIds());
                    size = imagesmd.size();
                    String question = String.format("<html>Are you sure you want to create a new set of %s with <font color=red>%s</font> %s?", type, size, (size > 1)?"elements":"element");
                    ScipionMessageDialog dlg = new ScipionMessageDialog(ScipionGalleryJFrame.this, "Question", question, msgfields);
                    int create = dlg.action;
                    if (create == ScipionMessageDialog.OK_OPTION) 
                    {
                        String[] command = new String[]{python, script, projectid, inputid, sqlitefile + "," + ((ScipionGalleryData)data).getPrefix(), String.format("SetOf%s", type), dlg.getFieldValue(runNameKey), other};
                        createSubset(command);
                    }
                }
            });
            if(data.hasClasses())
            {
                classcmdbutton = getScipionButton("Create Classes", new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        MetaData md = data.getMd(data.getEnabledIds());
                        int size = md.size();
                        String msg = String.format("<html>Are you sure you want to create a new set of Classes with <font color=red>%s</font> %s?", size, (size > 1)?"elements":"element");
                        ScipionMessageDialog dlg = new ScipionMessageDialog(ScipionGalleryJFrame.this, "Question", msg, msgfields);
                        int create = dlg.action;
                        if (create == ScipionMessageDialog.OK_OPTION) {
                            String output = ((ScipionGalleryData)data).getSelf().equals("Class2D")? "SetOfClasses2D":"SetOfClasses3D";
                            String[] command = new String[]{python, script, projectid, inputid, sqlitefile + "," + ((ScipionGalleryData)data).getPrefix(), output , dlg.getFieldValue(runNameKey), other};
                            createSubset(command);
                            
                        }

                    }
                });
                
                buttonspn.add(classcmdbutton);
            }
            
            buttonspn.add(cmdbutton);
            pack();
            enableActions();
            jcbBlocks.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    enableActions();
                }
            });
        }

    }

    public void reloadTableData(boolean changed)
    {
        super.reloadTableData(changed);
        enableActions();
    }
    

    public JButton getScipionButton(String text, ActionListener listener) {

        JButton button = new JButton(text);
        button.addActionListener(listener);

        return button;
    }

    

    protected void enableActions() {
        boolean isenabled = data.allowGallery();
        Color color = isenabled ? ScipionMessageDialog.firebrick : ScipionMessageDialog.lightgrey;
        Color forecolor = isenabled ? Color.WHITE : Color.GRAY;
        if(cmdbutton != null)
        {
            cmdbutton.setEnabled(isenabled);
            cmdbutton.setBackground(color);
            cmdbutton.setForeground(forecolor);
        }
        if(classcmdbutton != null)
        {
            isenabled = data.hasClasses() && !data.isVolumeMode();
            color = isenabled? ScipionMessageDialog.firebrick: ScipionMessageDialog.lightgrey; 
            forecolor = isenabled? Color.WHITE: Color.GRAY;
            classcmdbutton.setEnabled( isenabled);
            classcmdbutton.setBackground(color);
            classcmdbutton.setForeground(forecolor);
        }
    }

  
    	public boolean proceedWithChanges()
	{
            return true;//without asking for changes
	}
    
   protected void createSubset(final String[] command) 
    {
        XmippWindowUtil.blockGUI(ScipionGalleryJFrame.this, "Creating set ...");
        new Thread(new Runnable() {

            @Override
            public void run() {

                try {
                    ((ScipionGalleryData)data).overwrite(sqlitefile);
                    String output = XmippUtil.executeCommand(command);
                    XmippWindowUtil.releaseGUI(ScipionGalleryJFrame.this.getRootPane());
                    if (output != null && !output.isEmpty()) {
                        System.out.println(output);
                        XmippDialog.showInfo(ScipionGalleryJFrame.this, output);
                        
                    }

                } catch (Exception ex) {
                    throw new IllegalArgumentException(ex.getMessage());
                }

            }
        }).start();
    }
        


}
