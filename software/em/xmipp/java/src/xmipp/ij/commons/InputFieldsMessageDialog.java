/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.ij.commons;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Dictionary;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import xmipp.utils.XmippWindowUtil;

/**
 *
 * @author airen
 */
public class InputFieldsMessageDialog extends JDialog implements ActionListener {
    
    private String msg;
    private JButton cancelbt;
    private JButton okbt;
    public int action;
   
    public static final int OK_OPTION = 1;
    public static final int CANCEL_OPTION = 0;
    private Map<String, String> fields;
    Map<String, JTextField> fieldstfs;
   
    
    public InputFieldsMessageDialog(JFrame parent, String title, String msg, Map<String, String> fields)
    {
        super(parent, true);
        this.msg = msg;
        this.fields = fields;
        fieldstfs = new HashMap<String, JTextField>();
        initComponents(title);
        
    }
    
     public InputFieldsMessageDialog(JFrame parent, String title, String msg)
    {
        this(parent, title, msg, null);
    }
     
    protected void initComponents(String title)
    {
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setTitle(title);
        setLayout(new GridBagLayout());
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.insets = new Insets(5, 5, 5, 5);
        JLabel msgLb = new JLabel(msg);
        add(msgLb, XmippWindowUtil.getConstraints(constraints, 0, 0, 2));
        int index = 1;
        JTextField tf;
        if(fields != null)
            for(Map.Entry<String, String> e: fields.entrySet())
            {

                add(new JLabel(e.getKey()), XmippWindowUtil.getConstraints(constraints, 0, index));
                tf = new JTextField(e.getValue());
                tf.setColumns(20);
                fieldstfs.put(e.getKey(), tf);
                add(tf, XmippWindowUtil.getConstraints(constraints, 1, index));
                index ++;
            }
        JPanel buttonspn = new JPanel();
        cancelbt = new JButton("Cancel");
        cancelbt.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                action = CANCEL_OPTION;
                close();
            }
        });
        buttonspn.add(cancelbt);
        okbt = new JButton("Ok");
        okbt.setBackground(XmippWindowUtil.firebrick);
        okbt.setForeground(Color.WHITE);
        okbt.addActionListener(this);
        
        buttonspn.add(okbt);
        add(buttonspn, XmippWindowUtil.getConstraints(constraints, 0, index, GridBagConstraints.HORIZONTAL));
        pack();
        if (getParent() != null)
            XmippWindowUtil.centerWindows(this, getParent());
        setVisible(true);
        
    }
    
 
    
    
    public void close()
    {
        setVisible(false);
                dispose();
    }

    @Override
    public void actionPerformed(ActionEvent ae) {
        action = OK_OPTION;
        close();
    }
    
    public String getFieldValue(String field)
    {
        return fieldstfs.get(field).getText();
    }
    
    
  
}
