/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import java.awt.Color;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.Filename;
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
    private String scripts, script, ctfscript;
    private String projectid;
    private JButton cmdbutton;
    private String sqlitefile;
    private String outputdir;
    private JButton classcmdbutton;
    private String python;
    private String inputid;
    private HashMap<String, String> msgfields;
    private final String runNameKey = "Run name:";
    private String other;
    private JButton representativesbt;
    ScipionMessageDialog dlg;
    
   

    public ScipionGalleryJFrame(String filename, ScipionMetaData md, ScipionParams parameters) {
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
            type = ((ScipionGalleryData)data).getScipionType() + "s";
            python = parameters.python;
            scripts = parameters.scripts;
            script = parameters.scripts + File.separator + "pw_create_image_subset.py";
            ctfscript = parameters.scripts + File.separator + "pw_recalculate_ctf.py";
            projectid = parameters.projectid;
            inputid = parameters.inputid;
            outputdir = new File(data.getFileName()).getParent();
            sqlitefile = String.format("%s%sselection%s", outputdir, File.separator, data.getFileExtension());
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
        Icon icon = new ImageIcon(Toolkit.getDefaultToolkit().getImage(Filename.getXmippPath("resources" + File.separator
						+ "fa-times.png")));
        JButton closebt = new JButton("Close", icon);
        closebt.addActionListener(new ActionListener() {

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
                    int size = 0;
                    
                    if(data.hasClasses())
                    {
                        for(ScipionMetaData.EMObject emo: ((ScipionGalleryData)data).getEnabledObjects())
                            if(emo.childmd != null)
                                size += emo.childmd.getEnabledObjects().size();
                    }
                    else
                        size = data.getEnabledIds().size();
                    if (confirmCreate(type, size)) 
                    {
                        String[] command = new String[]{python, script, projectid, inputid, sqlitefile + "," + ((ScipionGalleryData)data).getPrefix(), String.format("SetOf%s", type), dlg.getFieldValue(runNameKey), other};
                        runCommand(command, "Creating set ...");
                    }
                }
            });
            if(data.hasClasses())
            {
                classcmdbutton = getScipionButton("Create Classes", new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        int size = data.getEnabledIds().size();
                       
                        if (confirmCreate("Classes", size)) {
                            String output = ((ScipionGalleryData)data).getSelf().equals("Class2D")? "SetOfClasses2D":"SetOfClasses3D";
                            String[] command = new String[]{python, script, projectid, inputid, sqlitefile + ",", output , dlg.getFieldValue(runNameKey), other};
                            runCommand(command, "Creating set ...");
                            
                        }

                    }
                });
                representativesbt = getScipionButton("Create Representatives", new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        MetaData md = data.getMd(data.getEnabledIds());
                        int size = md.size();
                        
                        if (confirmCreate("Representatives", size)) {
                            String output = ((ScipionGalleryData)data).getSelf().equals("Class2D")? "SetOfParticles,Representatives":"SetOfVolumes,Representatives";
                            String[] command = new String[]{python, script, projectid, inputid, sqlitefile + ",", output , dlg.getFieldValue(runNameKey), other};
                            runCommand(command, "Creating set ...");
                            
                        }

                    }
                });
                
                buttonspn.add(representativesbt);
                buttonspn.add(classcmdbutton);
                
            }
            
            buttonspn.add(cmdbutton);
            if(data.isCTFMd())
            {
                icon = new ImageIcon(Toolkit.getDefaultToolkit().getImage(Filename.getXmippPath("resources" + File.separator
						+ "fa-cogs.png")));
                JButton recalculatectfbt = getScipionButton("Recalculate CTFs", new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        String recalculatefile = outputdir + File.separator + "ctfrecalculate.txt";
                        ((ScipionGalleryData)data).exportCTFRecalculate(recalculatefile);
                        String[] command = new String[]{python, ctfscript, projectid, inputid, recalculatefile};
                        runCommand(command, "Recalculating CTFs ...");
                    }
                }, icon);
                buttonspn.add(recalculatectfbt);
            }
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
    
    public boolean confirmCreate(String output, int size)
    {
        String msg = String.format("<html>Are you sure you want to create a new set of %s with <font color=red>%s</font> %s?", output, size, (size > 1)?"elements":"element");
        if(data.size() == size)
            msg += "<br><font color=red>Note:</font> There are no disabled items to dismiss";
        dlg = new ScipionMessageDialog(ScipionGalleryJFrame.this, "Question", msg, msgfields);
                        int create = dlg.action;
        return (create == ScipionMessageDialog.OK_OPTION);
    }

    public void reloadTableData(boolean changed)
    {
        super.reloadTableData(changed);
        enableActions();
    }
    
    
    public JButton getScipionButton(String text, ActionListener listener) {
        Image imp = Toolkit.getDefaultToolkit().getImage(Filename.getXmippPath("resources" + File.separator
						+ "fa-plus-circle.png"));
        Icon icon = new ImageIcon(imp);
        return getScipionButton(text, listener, icon);
    }
    
    public JButton getScipionButton(String text, ActionListener listener, Icon icon) {
        JButton button = new JButton(text.replace("Create ", ""), icon);
        button.setToolTipText(text);
        button.addActionListener(listener);
        button.setBackground(ScipionMessageDialog.firebrick);
        button.setForeground(Color.WHITE);
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
            representativesbt.setEnabled( isenabled);
            representativesbt.setBackground(color);
            representativesbt.setForeground(forecolor);
        }
    }

  
    public boolean proceedWithChanges()
    {
        return true;//without asking for changes
    }
    
   protected void runCommand(final String[] command, String msg) 
    {
        XmippWindowUtil.blockGUI(ScipionGalleryJFrame.this, msg);
        new Thread(new Runnable() {

            @Override
            public void run() {

                try {
                    ((ScipionGalleryData)data).overwrite(sqlitefile);
                    String output = XmippUtil.executeCommand(command);
                    XmippWindowUtil.releaseGUI(ScipionGalleryJFrame.this.getRootPane());
//                    if (output != null && !output.isEmpty()) {
//                        System.out.println(output);
//                        XmippDialog.showInfo(ScipionGalleryJFrame.this, output);
//                        
//                    }

                } catch (Exception ex) {
                    ex.printStackTrace();
                    throw new IllegalArgumentException(ex.getMessage());
                }

            }
        }).start();
    }
        
   

    /**
	 * Open another metadata separataly *
	 */
    @Override
	public void openMetadata(MetaData md)
	{
            try
            {
                String[] args = new String[]{"--scipion", python, scripts, projectid, inputid, ""};
                ScipionParams params = new ScipionParams(args);
                new ScipionGalleryJFrame(new ScipionGalleryData(this, params, (ScipionMetaData)md));
            }
            catch(Exception e)
            {
                XmippDialog.showError(this, e.getMessage());
            }
	}

}
