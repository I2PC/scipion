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
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import xmipp.jni.Filename;
import xmipp.jni.MetaData;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippWindowUtil;

/**
 *
 * @author airen
 */
public class ExportImagesJDialog extends JDialog{
    private MetaData md;
    private int label;
    private String path;
    private XmippFileChooser fc;
    private JButton savebt;
    private JButton cancelbt;
    private String note = "<html> <font color='red'>Extensions supported are: jpg, mrc, ...</font> ";
    private JButton browsebt;
    private JTextField pathtf;
    
    public ExportImagesJDialog(GalleryJFrame parent)
    {
        super(parent);
        
        md = parent.data.md;
        path = parent.data.getFileName();
        path = Filename.removeExtension(path) + "_export.stk";
        label = parent.data.getRenderLabel();
        initComponents();
    }

    private void initComponents() {
        setTitle("Export Images ...");
        
        setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.insets = new Insets(5, 5, 5, 5);
        add(new JLabel("Path:"), XmippWindowUtil.getConstraints(c, 0, 0));
        

        pathtf = new JTextField(path);
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
        add(new JLabel(note), XmippWindowUtil.getConstraints(c, 0, 1, GridBagConstraints.HORIZONTAL));
        
        JPanel actionspn = new JPanel();
        cancelbt = XmippWindowUtil.getTextButton("Cancel", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                close();
            }
        });
        actionspn.add(cancelbt);
                
        savebt = XmippWindowUtil.getTextButton("Save", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                saveImages();
                close();
            }
        });
        actionspn.add(savebt);
        add(actionspn, XmippWindowUtil.getConstraints(c, 0, 2, GridBagConstraints.HORIZONTAL));
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
        path = pathtf.getText();
        md.writeImages(path, false, label);
    }
    
}
