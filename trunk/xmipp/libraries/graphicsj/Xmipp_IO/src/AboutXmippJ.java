
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
public class AboutXmippJ extends PlugInFrame {

    public AboutXmippJ() {
        super("About XmippJ...");

        JLabel jlIcon = new JLabel(new ImageIcon(getClass().getResource("/resources/i2pc.png")));
        JTextPane jtpAbout = new JTextPane();

        jtpAbout.setContentType("text/html");
        jtpAbout.setEditable(false);
        jtpAbout.setText("<html>"
                + "Juanjo Vega (<b>juanjo.vega@gmail.com</b>)<br>"
                + "<hr>"
                + "Instruct Image Processing Center.<br>"
                + "Biocomputing Unit.<br>"
                + "National Center for Biotechnology (CNB/CSIC).<br>"
                + "2010</html>");

        add(jlIcon, java.awt.BorderLayout.WEST);
        add(jtpAbout, java.awt.BorderLayout.CENTER);

        pack();
        setLocationRelativeTo(null);
        setResizable(false);
    }

    @Override
    public void run(String arg) {
        setVisible(true);
    }
}
