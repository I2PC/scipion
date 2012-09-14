package xmipp.ij.plugins.maskstoolbar;

import ij.plugin.frame.PlugInFrame;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class AboutMasksToolBar extends PlugInFrame {

    public AboutMasksToolBar() {
        super("About Masks Tool Bar...");

        JLabel jlIcon = new JLabel(new ImageIcon(getClass().getResource("/resources/I2PC.png")));
        JScrollPane jspText = new JScrollPane();
        JTextPane jtpAbout = new JTextPane();

        jtpAbout.setContentType("text/html");
        jtpAbout.setEditable(false);
        jtpAbout.setText("<html>\nJuanjo Vega (<b>juanjo.vega@gmail.com</b>)<br>"
                + "Biocomputing Unit.<br>"
                + "National Center for Biotechnology (CNB/CSIC).<br>"
                + "Madrid. May 2011.</html>");
        jspText.setViewportView(jtpAbout);

        add(jlIcon, java.awt.BorderLayout.WEST);
        add(jspText, java.awt.BorderLayout.CENTER);

        pack();
        setLocationRelativeTo(null);
        setResizable(false);
    }

    @Override
    public void run(String arg) {
        setVisible(true);
    }
}
