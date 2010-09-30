
import browser.ICONS_MANAGER;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.frame.PlugInFrame;
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
public class Main extends PlugInFrame {

    JLabel label;

    public Main() {
        super("[Xmipp IO - TEST]");

        label = new JLabel("");
        add(label);

        setLocationRelativeTo(null);
        setVisible(true);
    }

    public static void main(String args[]) {
        new ImageJ();

        (new Main()).run("");
    }

    @Override
    public void run(String arg) {
        int W = ICONS_MANAGER.PREVIEW_WIDTH;
        int H = ICONS_MANAGER.PREVIEW_HEIGHT;
        String dir = "//home//juanjo//temp//";
        String file = "Iter_6_filtered_reconstruction.vol";

        Spider_Reader sr = new Spider_Reader();

        ImagePlus ip = sr.loadThumbnail(dir, file, W, H);

        label.setIcon(new ImageIcon(ip.getImage()));

        pack();

        super.run(arg);
    }
}
