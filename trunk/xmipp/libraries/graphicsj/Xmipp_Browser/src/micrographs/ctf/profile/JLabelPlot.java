package micrographs.ctf.profile;


import ij.ImagePlus;
import javax.swing.ImageIcon;
import javax.swing.JLabel;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class JLabelPlot extends JLabel {

    ImageIcon icon = new ImageIcon();

    public JLabelPlot() {
        super();
        icon = new ImageIcon();

        setIcon(icon);
    }

    public void setImage(ImagePlus ip) {
        icon.setImage(ip.getImage());
        repaint();
    }
}
