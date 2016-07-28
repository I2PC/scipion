/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import xmipp.ij.commons.InputFieldsMessageDialog;
import xmipp.utils.ScipionParams;
import java.awt.Color;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.SwingUtilities;
import xmipp.ij.commons.XmippApplication;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryJFrame extends GalleryJFrame {

    private static Logger logger = Logger.getLogger(ScipionGalleryJFrame.class.getName());

    private String type;
    private Integer port;
    private JButton cmdbutton;
    private String sqlitefile;
    private JButton classcmdbutton;
    private String inputid;
    private HashMap<String, String> msgfields;
    private final String runNameKey = "Run name:";
    private String other;
    private JButton representativesbt;
    private InputFieldsMessageDialog dlg;
    private JButton createvolbt;
    private String setType;
    private boolean recalculateCTF;

    private static final String runProtCreateSubset = "run protocol ProtUserSubSet inputObject=%s sqliteFile='%s','%s' outputClassName=%s other='%s' label='%s'";
    
   
    public ScipionGalleryJFrame(ScipionGalleryData data) {
        super(data);
        readScipionParams((ScipionParams)data.parameters);
        
    }
      
  

    protected void readScipionParams(ScipionParams parameters)
    {
        try {
            setType = ((ScipionMetaData)data.getMd()).getSetType();
            type = data.hasClasses()? "Particles": setType.replace("SetOf", "");
            port = parameters.port;
            inputid = parameters.inputid;
            sqlitefile = data.getTmpFile("_state");
            msgfields = new HashMap<String, String>();
            msgfields.put(runNameKey, "create subset");
            other = parameters.other;
            recalculateCTF = parameters.recalculateCTF;
            initComponents();
        } catch (Exception ex) {
            logger.log(Level.SEVERE, null, ex);
            throw new IllegalArgumentException(ex.getMessage());
        }
    }
    
    protected void initComponents() {
        Icon icon = XmippResource.getIcon("fa-times.png");
        JButton closebt = new JButton("Close", icon);
        closebt.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                close();
            }
        });
        
        buttonspn.add(closebt);
        if(!XmippApplication.isScipion())
            return;
            
        if (type != null) {
            if(!data.isCTFMd())
            {
                cmdbutton = XmippWindowUtil.getScipionIconButton("Create " + type);
                cmdbutton.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createSimpleSubset();
                    }
                });
                buttonspn.add(cmdbutton);
            }
            if(data.hasClasses())
            {
                classcmdbutton = XmippWindowUtil.getScipionIconButton("Create Classes");
                classcmdbutton.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createSubsetFromClasses();

                    }
                });
                
                String repText = setType.equals("SetOfClasses2D") ? "Create Averages": "Create Volumes";
                representativesbt = XmippWindowUtil.getScipionIconButton(repText);
                representativesbt.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createRepresentativesSubset();

                    }
                });
                
                buttonspn.add(representativesbt);
                buttonspn.add(classcmdbutton);
            }
            
            
            if(data.isCTFMd())
            {
                if (recalculateCTF) {
                    icon = XmippResource.getIcon("fa-cogs.png");
                    JButton recalculatectfbt = XmippWindowUtil.getScipionIconButton("Recalculate CTFs");
                    recalculatectfbt.addActionListener(new ActionListener() {

                        @Override
                        public void actionPerformed(ActionEvent ae) {
                            if (!data.hasRecalculateCTF()) {
                                XmippDialog.showError(ScipionGalleryJFrame.this, "There are no ctfs to recalculate");
                                return;
                            }
                            String command = String.format("run function recalculateCTF %s %s", inputid, sqlitefile);
                            runCommand(command, "Recalculating CTF", false);
                        }
                    });
                    recalculatectfbt.setIcon(icon);
                    buttonspn.add(recalculatectfbt);
                }

                JButton ctfsubsetbt = XmippWindowUtil.getScipionIconButton("Create Micrographs");
                ctfsubsetbt.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        if (confirmCreate("Micrographs")) 
                        {
                            String command = String.format(runProtCreateSubset, 
                            inputid, sqlitefile, "", "SetOfMicrographs", other, getRunLabel());
                            runCommand(command, "Creating set ...");
                        }
                    }
                });
                buttonspn.add(ctfsubsetbt);

            }
            if(setType.equals("SetOfVolumes") || setType.equals("SetOfClasses3D"))
            {
                createvolbt = XmippWindowUtil.getScipionIconButton("Create Volume");
                createvolbt.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createVolume();
                    }
                });
                buttonspn.add(createvolbt);
                createvolbt.setVisible(!data.isTableMode());
            }


            if (!SwingUtilities.isEventDispatchThread()) {

                Runnable pack = new Runnable() {
                    @Override
                    public void run() {
                        ScipionGalleryJFrame.this.pack();
                    }
                };

                try {
                    SwingUtilities.invokeAndWait(pack);
                } catch (Exception e) {

                    logger.log(Level.WARNING, "ScipionGalleryJFrame.pack threw an exception: " + e.getMessage());
                }
            } else {
                pack();
            }

            enableActions();
            jcbBlocks.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    enableActions();
                }
            });
        }
    }
    
    
    protected void createVolume()
    {
        msgfields.put(runNameKey, "ProtRegisterVolume");
        dlg = new InputFieldsMessageDialog(ScipionGalleryJFrame.this, "Question", "Are you sure you want to register volume in scipion?", msgfields);
        int create = dlg.action;
        if (create == InputFieldsMessageDialog.OK_OPTION);
        {
            String command = String.format(runProtCreateSubset, 
                inputid, sqlitefile, "", setType, data.getSelVolId().toString() + ",Volume", getRunLabel());
            runCommand(command, "Creating volume ...");
        }
    }
    
    public String getRunLabel()
    {
        return dlg.getFieldValue(runNameKey);
    }
    
    public String getDBPreffix()
    {
        return ((ScipionGalleryData)data).getPreffix();
    }
    
    protected void createSimpleSubset()
    {
        int size = 0;
       
        if(data.hasClasses())
        {
            boolean[] selection = gallery.getSelection();
            for(ScipionMetaData.EMObject emo: ((ScipionGalleryData)data).getEMObjects())
            {
                if(gallery.hasSelection() && !data.hasDisabled())
                { 
                	if(selection[emo.index] && emo.childmd != null)
                		size += emo.childmd.getEnabledCount();
                }
                else if(emo.isEnabled() && emo.childmd != null)
                    size += emo.childmd.getEnabledCount();
            }
        }
        else if (gallery.hasSelection() && !data.hasDisabled())
            size = gallery.getSelectionCount();
        else
            size = data.getEnabledCount();
        
        if (confirmCreate(type, size)) 
        {
            String command = String.format(runProtCreateSubset, 
                    inputid, sqlitefile, data.getPreffix(), String.format("SetOf%s", type), other, getRunLabel());
            runCommand(command, "Creating set ...");
        }
    }

    protected void createSubsetFromClasses()
    {
        if (confirmCreate("Classes")) {
            String command = String.format(runProtCreateSubset, 
                inputid, sqlitefile, "", setType , other, getRunLabel());
            runCommand(command, "Creating set ...");
        }
    }
    
    protected void createRepresentativesSubset()
    {
        if (confirmCreate("Representatives")) {
            String output = setType.equals("SetOfClasses2D")? "SetOfAverages,Representatives":"SetOfVolumes,Representatives";
            String command = String.format(runProtCreateSubset, inputid, sqlitefile, "", output , other, getRunLabel());
            runCommand(command, "Creating set ...");
        }
    }
    
    public boolean confirmCreate(String output)
    {
    	int size;
    	if(!data.hasDisabled())
    	{
    		size = gallery.hasSelection()? gallery.getSelectionCount(): data.size();
    	}
    	else
    		size = data.getEnabledCount();
        return confirmCreate(output, size);
    }
    
    public boolean confirmCreate(String output, int size)
    {
        String msg = String.format("<html>Are you sure you want to create a new set of %s with <font color=red>%s</font> %s?", output, size, (size == 1)?"element":"elements");
        dlg = new InputFieldsMessageDialog(ScipionGalleryJFrame.this, "Question", msg, msgfields);
        return dlg.action == InputFieldsMessageDialog.OK_OPTION;
    }

    public void reloadTableData(boolean changed)
    {
        super.reloadTableData(changed, gallery.getSelection());
        enableActions();
    }

    protected void enableActions() {
        boolean isenabled = !data.isVolumeMode();
        Color color = isenabled ? XmippWindowUtil.firebrick : XmippWindowUtil.lightgrey;
        Color forecolor = isenabled ? Color.WHITE : Color.GRAY;
        if(cmdbutton != null)
        {
            cmdbutton.setVisible(isenabled);
            cmdbutton.setBackground(color);
            cmdbutton.setForeground(forecolor);
        }
        if(classcmdbutton != null)
        {
            isenabled = data.hasClasses() && !data.isVolumeMode();
            color = isenabled? XmippWindowUtil.firebrick: XmippWindowUtil.lightgrey; 
            forecolor = isenabled? Color.WHITE: Color.GRAY;
            classcmdbutton.setVisible( isenabled);
            classcmdbutton.setBackground(color);
            classcmdbutton.setForeground(forecolor);
            representativesbt.setVisible( isenabled);
            representativesbt.setBackground(color);
            representativesbt.setForeground(forecolor);
        }
    }

    @Override
    protected void changeView()
    {
        super.changeView();
        
        if(setType.equals("SetOfVolumes") || setType.equals("SetOfClasses3D"))
        {

            // It will be null when invoked as a stand alone viewer.
        	if (cmdbutton != null){
        		cmdbutton.setVisible(data.isTableMode());
            	createvolbt.setVisible(!data.isTableMode());
        	}
        }
    }
  
    public boolean proceedWithChanges()
    {
        return true;
    }
    
    protected void runCommand(final String command, String msg) 
    {
        runCommand(command, msg, gallery.hasSelection() && !data.hasDisabled());
    }
    
    protected void runCommand(final String command, String msg, boolean useSelection) 
    {
        try {
            boolean[] selection = null;
            
            if(useSelection)
                selection = gallery.getSelection();
            ((ScipionGalleryData)data).overwrite(sqlitefile, selection);
            XmippWindowUtil.runCommand(command, port);
            //close(false);
        } catch (SQLException ex) {
            logger.log(Level.SEVERE, null, ex);
        }
    }
    /**
	 * Open another metadata separataly *
	 */
    @Override
    public void openMetadata(final MetaData md)
    {
        try
        {
            SwingUtilities.invokeLater(new Runnable() {

                @Override
                public void run() {
                    new ScipionGalleryJFrame(new ScipionGalleryData(ScipionGalleryJFrame.this, ((ScipionParams)data.parameters).getScipionParams(), (ScipionMetaData)md));
                }
            });
            
        }
        catch(Exception e)
        {
            XmippDialog.showError(this, e.getMessage());
        }
    }
    
    protected void initGalleryMenu() {
            menu = new ScipionGalleryMenu();
                    
    }
    
    protected class ScipionGalleryMenu extends GalleryMenu//To customize showj menu for scipion
    {
        
        @Override
        protected void handleActionPerformed(ActionEvent evt)
        {
            super.handleActionPerformed(evt);
            String cmd = evt.getActionCommand();
            try
            {
                    if (cmd.equals(FILE_LOAD_SEL))
                    {
                    	fc.setApproveButtonText("Open");
                        fc.setDialogTitle("Load a status file");
                        fc.setApproveButtonToolTipText("Choose a file with state data");

                        if (fc.showOpenDialog(ScipionGalleryJFrame.this) != XmippFileChooser.CANCEL_OPTION)
                            loadSelection(fc.getSelectedPath());
                    }
                    if (cmd.equals(FILE_SAVE_SEL))
                    {
                        fc.setSelectedFile(new File(sqlitefile));
                        fc.setDialogTitle("Save state");
                        fc.setApproveButtonToolTipText("Choose where to save the status");
                        fc.setApproveButtonText("Save");
                         if (fc.showOpenDialog(ScipionGalleryJFrame.this) != XmippFileChooser.CANCEL_OPTION)
                            saveSelection(fc.getSelectedPath());
                    }
            }
            catch (Exception e)
            {
                    showException(e);
            }
        }

        protected void loadSelection(String path) {
            
                ((ScipionGalleryData)data).loadSelection(path);
                reloadTableData();
            
        }

        protected void saveSelection(String path) {
            try {
//                boolean[] selection = null;
//                if(gallery.hasSelection() && !data.hasDisabled())
//                    selection = gallery.getSelection();
//                ((ScipionGalleryData)data).overwrite(path, selection);

                ScipionGalleryData sgData = (ScipionGalleryData)data;

                if (sgData.getScipionMetaData().getParent() != null) {
                    sgData.getScipionMetaData().getParent().overwrite(getFilename(), path, null);
                } else {
                    sgData.overwrite(path, null);
                }
            } catch (SQLException ex) {
                logger.log(Level.SEVERE, null, ex);
            }
        }

		@Override
		public void addExtraMenuItems()
		{
			addItem(FILE_LOAD_SEL, "Load state ...");
            addItem(FILE_SAVE_SEL, "Save state ...", "save_as.gif");
			
		}
    }
        
        

}
