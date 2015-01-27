/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.process.ByteProcessor;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JToggleButton;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.ScipionParams;
import xmipp.viewer.scipion.ScipionGalleryJFrame;
import xmipp.viewer.scipion.ScipionMessageDialog;

/**
 *
 * @author airen
 */
public class MaskJFrame extends JFrame{
    protected ImagePlus mask;
    private ImageJPanel imagepn;
    private JButton registerbt;
    private JButton previewbt;
    protected int width, height;
    private ImagePlus imp;
    private JToggleButton invertbt;
    private JToggleButton addbt;
    private final XmippImageWindow iw;
    
    public MaskJFrame(XmippImageWindow iw)
    {
        this.iw = iw;
        createMask();
        XmippUtil.showImageJ(Tool.VIEWER);
        IJ.setTool(Toolbar.OVAL);
        initComponents();
    }
    
    protected void initComponents()
    {
        setTitle("Mask Manager");
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
       
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.insets = new Insets(0, 5, 0, 5);
        constraints.anchor = GridBagConstraints.WEST;
        setLayout(new GridBagLayout());
        
        imagepn = new ImageJPanel(mask, width, height);
        add(imagepn, XmippWindowUtil.getConstraints(constraints, 0, 0));
        invertbt = new JCheckBox("Invert");
        invertbt.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                createMask();
                imagepn.setMask(mask, width, height);
            }

            
        });
        add(invertbt);
        add(invertbt, XmippWindowUtil.getConstraints(constraints, 0, 1));
        JPanel commandspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        add(commandspn, XmippWindowUtil.getConstraints(constraints, 0, 2, 1, 1, GridBagConstraints.HORIZONTAL));
        
              
        registerbt = XmippWindowUtil.getScipionButton("Create Mask");
        registerbt.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                registerMask();
            }
            
        });
        
        commandspn.add(registerbt);
        registerbt.setVisible(XmippApplication.isScipion());
        pack();
        setVisible(true);

    }
    
    protected void registerMask() {
        HashMap<String, String> msgfields = new HashMap<String, String>();
        String field = "Run name:";
        msgfields.put(field, "create mask");
        ScipionMessageDialog dlg = new ScipionMessageDialog(this, "Question", "Are you sure you want to register mask?", msgfields);
                        int create = dlg.action;
        boolean register = (create == ScipionMessageDialog.OK_OPTION);
        if (register)
            try {
            String path = "mask.spi";
            ImageGeneric ig = XmippImageConverter.convertToImageGeneric(mask);
            ig.write(path);
            ScipionParams params = (ScipionParams)iw.getParams();
            String label = dlg.getFieldValue(field);
            String[] command = new String[]{params.python, params.getRegisterMaskScript(), params.projectid, params.inputid, path, label};

            String output = XmippWindowUtil.executeCommand(command, true);
            System.out.println(output);
            } catch (Exception ex) {
                ex.printStackTrace();
                Logger.getLogger(MaskJFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
    }
    
   
    
    private void createMask() {
            imp = IJ.getImage();
            width = imp.getCanvas().getWidth();
            height = imp.getCanvas().getHeight();

            if (imp != null) {
                Roi roi = imp.getRoi();

                if (roi != null) {
                    if (mask == null || !mask.isVisible()) {
                        createDefaultMask();
                    }

                    mask.getProcessor().setColor(Color.WHITE);
                    mask.getProcessor().fill(roi);
                    mask.updateAndDraw();
                    if(isInvert())
                    {
                        mask.setRoi((Roi) null); // ...clears selection...
                        mask.getProcessor().invert();
                        mask.setRoi(roi);   // ...to restore it later.
                        mask.updateAndDraw();
                    }
                    
                } else 
                    createDefaultMask();
                    
                
            } else {
                IJ.error("There are no images open.");
            }
            
    }
    
    public boolean isInvert()
    {
        if(invertbt == null)
            return false;
        return invertbt.isSelected();
    }
    
    public boolean isAddRoi()
    {
        if(addbt == null)
            return false;
        return addbt.isSelected();
    }
    
    protected void createDefaultMask()
    {
        ByteProcessor processor = new ByteProcessor(imp.getWidth(), imp.getHeight());
        mask = new ImagePlus("Mask", processor);
        mask.updateAndDraw();
            

    }
    
    public void refreshMask()
    {
        createMask();
        imagepn.setMask(mask, width, height);
        pack();
    }
    
        
   
   
    
}
