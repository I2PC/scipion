/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.windows;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import xmipp.ij.commons.XmippApplication;
import xmipp.ij.commons.XmippUtil;
import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippWindowUtil;

/**
 *
 * @author airen
 */
public class ExportImagesJDialog extends JDialog{

    private int label;
    private String path;
    private XmippFileChooser fc;
    private JButton savebt;
    private JButton cancelbt;
    private final String note1 = "<html><b>Note 1:</b> Only enabled images will be saved";
    private final String note2 = "<html><b>Note 2:</b> Use extension stk, mrcs or img to save as SPIDER, MRC or IMAGIC stacks";
    private JButton browsebt;
    private JTextField pathtf;
    private JCheckBox applygeochb;
 
    private GalleryJFrame frame;
    
    public ExportImagesJDialog(GalleryJFrame parent)
    {
        super(parent);
        this.frame = parent;

        path = parent.data.getFileName();
        if(path == null)
            path = "images";
        path = Filename.removeExtension(path) + "_export.stk";
        label = parent.data.getRenderLabel();
        initComponents();
    }

    private void initComponents() {
        setTitle("Export Images ...");
        
        setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.anchor = GridBagConstraints.WEST;
        c.insets = new Insets(5, 5, 5, 5);
        add(new JLabel("Path:"), XmippWindowUtil.getConstraints(c, 0, 0));
        

        pathtf = new JTextField(path);
        pathtf.setColumns(50);
        add(pathtf, XmippWindowUtil.getConstraints(c, 1, 0));
        fc = new XmippFileChooser(path);
        if(path.contains(File.separator))
        {
            int index = path.lastIndexOf(File.separator);//-1 if separator does not exists
            String dir = path.substring(0, index);
            fc.setCurrentDirectory(new File(dir));
        }
        browsebt = XmippWindowUtil.getIconButton("folderopen.gif", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                int returnVal = fc.showOpenDialog(null);

		if (returnVal == XmippFileChooser.APPROVE_OPTION)
		{
			File file = fc.getSelectedFile();
			String text = file.getPath(); 
			pathtf.setText(text);
		}
            }
        });
        add(browsebt, XmippWindowUtil.getConstraints(c, 2, 0));
        JLabel applygeolb = new JLabel("Apply Geometry:"); 
        add(applygeolb, XmippWindowUtil.getConstraints(c, 0, 1));
        

        applygeochb = new JCheckBox();
        applygeochb.setSelected(frame.data.useGeo());
        applygeochb.setEnabled(frame.data.containsGeometryInfo());
        applygeolb.setEnabled(frame.data.containsGeometryInfo());
        add(applygeochb, XmippWindowUtil.getConstraints(c, 1, 1));
        add(new JLabel(note1), XmippWindowUtil.getConstraints(c, 0, 2, GridBagConstraints.HORIZONTAL));
        add(new JLabel(note2), XmippWindowUtil.getConstraints(c, 0, 3, GridBagConstraints.HORIZONTAL));
        JPanel actionspn = new JPanel();
        cancelbt = XmippWindowUtil.getTextButton("Cancel", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                close();
            }
        });
        actionspn.add(cancelbt);
                
        savebt = XmippWindowUtil.getTextButton("Export", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                saveImages();
                close();
            }
        });
        actionspn.add(savebt);
        add(actionspn, XmippWindowUtil.getConstraints(c, 0, 4, GridBagConstraints.HORIZONTAL));
        pack();
        XmippWindowUtil.setLocation(0.5, 0.5, this);
        setVisible(true);
    }
    
    private void close()
    {
        setVisible(false);
        dispose();
    }
    
    private void saveImages()
    {
            XmippWindowUtil.blockGUI(frame, "Exporting images..", true);
            new Thread(new Runnable() {

                @Override
                public void run() {

                    try {
                    	path = pathtf.getText();
                    	if(path.endsWith(".mrcs"))
                    		path += ":mrcs";
                    	int label = frame.data.getRenderColumn().label;
                    	MetaData md;
                    	if(frame.data.isScipionInstance())
                    	{
                    		md = frame.data.getImagesMd();//reads only enabled objects
                    		label = MDLabel.MDL_IMAGE;
                    	}
                    	else
                    		md = frame.data.getMd();
                    	
                    	md.writeMdToStack(path, applygeochb.isSelected(), frame.data.isWrap(), label);

                    } catch (Exception ex) {

                        XmippDialog.showException(frame,ex);
                    } finally {
                        XmippWindowUtil.releaseGUI(frame.getRootPane());
                    }

                }
            }).start();
    }
    
    

    
}
