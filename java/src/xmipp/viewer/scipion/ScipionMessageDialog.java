/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.scipion;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import xmipp.utils.XmippWindowUtil;

/**
 *
 * @author airen
 */
public class ScipionMessageDialog extends JDialog implements ActionListener {
    private String msg;
    private JFrame frame;
    private JButton cancelbt;
    private JButton okbt;
    private int action;
    public static final String firebrick = "#B22222";
    public static final String lightgrey = "#EAEBEC";
    public static final int OK_OPTION = 1;
    public static final int CANCEL_OPTION = 0;
    
    
    public ScipionMessageDialog(JFrame parent, String title, String msg)
    {
        super(parent, true);
        this.msg = msg;
        this.frame = parent;
        initComponents(title);
        
    }
    
    protected void initComponents(String title)
    {
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setTitle(title);
        setLayout(new GridBagLayout());
        GridBagConstraints constraints = new GridBagConstraints();
        JLabel msgLb = new JLabel(msg);
        add(msgLb, XmippWindowUtil.getConstraints(constraints, 0, 0));
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
        okbt = getScipionButton("Ok", true);
        okbt.addActionListener(this);
        
        buttonspn.add(okbt);
        add(buttonspn, XmippWindowUtil.getConstraints(constraints, 0, 1));
        pack();
        XmippWindowUtil.setLocation(0.5, 0.5, this);
        setVisible(true);
        
    }
    
    public static JButton getScipionButton(String text, boolean isenabled)
    {
        JButton button = new JButton(text);
        Color color = Color.decode(isenabled? ScipionMessageDialog.firebrick: ScipionMessageDialog.lightgrey); 
        Color forecolor = isenabled? Color.WHITE: Color.GRAY;
        button.setBackground(color);
        button.setForeground(forecolor);
        return button;
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
    
    
    public static int showQuestion(JFrame parent, String msg)
    {
        ScipionMessageDialog dlg = new ScipionMessageDialog(parent, "Question", msg);
        return dlg.action;
    }
    
}
