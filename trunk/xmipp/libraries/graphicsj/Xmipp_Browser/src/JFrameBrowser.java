
import browser.JPanelBrowser;
import ij.gui.GUI;
import java.awt.BorderLayout;
import javax.swing.JFrame;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class JFrameBrowser extends JFrame {

    public JFrameBrowser(String title) {
        super(title);
        setLayout(new BorderLayout());
        //setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JPanelBrowser jPanelBrowser = new JPanelBrowser();
        add(jPanelBrowser, BorderLayout.CENTER);

        pack();
        GUI.center(this);
//        setLocationRelativeTo(null);
//        setVisible(true);
    }
}
