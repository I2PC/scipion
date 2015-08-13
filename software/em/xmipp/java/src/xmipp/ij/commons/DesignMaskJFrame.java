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
import ij.process.FloatProcessor;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;

import xmipp.jni.ImageGeneric;
import xmipp.utils.ScipionParams;
import xmipp.utils.XmippWindowUtil;

/**
 *
 * @author airen
 */
public class DesignMaskJFrame extends JFrame implements ActionListener{
    protected ImagePlus mask;
    private ImageJPanel imagepn;
    private JButton registerbt;
    protected int width, height;
    private ImagePlus imp;
    private JCheckBox invertbt;
    private JCheckBox addbt;
    private final XmippImageWindow iw;
    private JToggleButton smoothbt;
    private JFormattedTextField smoothtf;
	private JButton refreshbt;
    
    public DesignMaskJFrame(XmippImageWindow iw)
    {
        this.iw = iw;
        createMask();
        XmippUtil.showImageJ(Tool.VIEWER);
        IJ.setTool(Toolbar.OVAL);
        initComponents();
    }
    
    protected void initComponents()
    {
        setTitle("Design Mask");
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
       
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.insets = new Insets(0, 5, 0, 5);
        constraints.anchor = GridBagConstraints.WEST;
        setLayout(new GridBagLayout());
        JToolBar optionspn = new JToolBar();

		optionspn.setFloatable(false);
        add(optionspn, XmippWindowUtil.getConstraints(constraints, 0, 0));
        refreshbt = XmippWindowUtil.getIconButton("fa-refresh.png", new ActionListener()
        {
        	
        	@Override
        	public void actionPerformed(ActionEvent e)
        	{
        		refreshMask();
        		
        	}
        });
        optionspn.add(refreshbt);
        invertbt = new JCheckBox("Invert:");
        invertbt.addActionListener(this);
        optionspn.add(invertbt);
        
        smoothbt = new JCheckBox("Smooth:");
        smoothbt.addActionListener(this);
        optionspn.add(smoothbt);
        smoothtf = new JFormattedTextField(NumberFormat.getIntegerInstance());
        smoothtf.setColumns(2);
        smoothtf.setValue(2);
        smoothtf.addActionListener(this);
        smoothtf.setEnabled(false);
        optionspn.add(smoothtf);
        imagepn = new ImageJPanel(mask, width, height);
        add(imagepn, XmippWindowUtil.getConstraints(constraints, 0, 1));
        
        
        JPanel commandspn = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        add(commandspn, XmippWindowUtil.getConstraints(constraints, 0, 2, 1, 1, GridBagConstraints.HORIZONTAL));
        
              
        registerbt = XmippWindowUtil.getScipionIconButton("Create Mask");
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
        InputFieldsMessageDialog dlg = new InputFieldsMessageDialog(this, "Question", "Are you sure you want to register mask?", msgfields);
                        int create = dlg.action;
        boolean register = (create == InputFieldsMessageDialog.OK_OPTION);
        if (register)
            try {
            String path = "mask.spi";
            ImageGeneric ig = XmippImageConverter.convertToImageGeneric(mask);
            ig.write(path);
            ScipionParams params = (ScipionParams)iw.getParams();
            String label = dlg.getFieldValue(field);
            String command = String.format("run protocol ProtCreateMask inputObj=%s maskFile='%s' label='%s'", params.inputid, path, label);
            XmippWindowUtil.runCommand(command, params.port);
//            String[] command = new String[]{params.python, params.getRegisterMaskScript(), params.projectid, params.inputid, path, label};
//            String output = XmippWindowUtil.executeCommand(command, true);
//            System.out.println(output);
            } catch (Exception ex) {
                ex.printStackTrace();
                Logger.getLogger(DesignMaskJFrame.class.getName()).log(Level.SEVERE, null, ex);
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
                if(isInvert())
                    mask.getProcessor().invert();
                if(isSmooth())
                    mask.getProcessor().blurGaussian(((Number) smoothtf.getValue()).intValue());
                
            } else 
                createDefaultMask();
           
            mask.getProcessor().multiply(1/mask.getProcessor().getMax());
            mask.getProcessor().setMinAndMax(0, 1);
            mask.updateAndDraw();    
            
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
    
    public boolean isSmooth()
    {
        if(smoothbt == null)
            return false;
        return smoothbt.isSelected();
    }
    
    public boolean isAddRoi()
    {
        if(addbt == null)
            return false;
        return addbt.isSelected();
    }
    
    protected void createDefaultMask()
    {
        FloatProcessor processor = new FloatProcessor(imp.getWidth(), imp.getHeight());
        mask = new ImagePlus("Mask", processor);
            

    }
    
    public void refreshMask()
    {
        createMask();
        imagepn.setMask(mask, width, height);
        pack();
    }

    @Override
    public void actionPerformed(ActionEvent ae) {
        smoothtf.setEnabled(isSmooth());
        createMask();
        imagepn.setMask(mask, width, height);
        
    }
    
        
   
   
    
}
