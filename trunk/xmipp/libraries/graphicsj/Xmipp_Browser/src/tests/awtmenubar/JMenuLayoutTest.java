package tests.awtmenubar;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JFrame;

/* Main-class to test the "dynamic" menu */
public class JMenuLayoutTest {

    public static void main(String args[]) {
        JFrame f = new JFrame("Test Layout");

        String[] menuLabels = {"File", "Edit", "Format", "View", "Tools", "Configure", "Admin", "Compile", "Options", "Help"};
        DynamicJMenuBar mbar = new DynamicJMenuBar(f, menuLabels);

        // make the menubar stand out
        //mbar.setBackground(Color.red);

        f.setMenuBar(mbar);
        f.setSize(800, 200);
        f.addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });
        f.setVisible(true);
    }
}
