

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import browser.DEBUG;
import browser.commandline.Parameters;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;
import micrographs.JFrameMicrographs;

/**
 *
 * @author Juanjo Vega
 */
public class XmippMicrographViewer implements PlugIn {

    public void run(String args) {
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "args" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            Parameters parameters = new Parameters(Macro.getOptions().split(" "));

            if (parameters.files != null) {
                for (int i = 0; i < parameters.files.length; i++) {
                    if (!parameters.files[i].isEmpty()) {
                        final String fileName = parameters.files[i];
                        java.awt.EventQueue.invokeLater(new Runnable() {

                            public void run() {
                                openFile(fileName);
                            }
                        });
                    }
                }
            }
        }
    }

    public static void openFile(String fileName) {
        try {
            JFrameMicrographs frameMicrographs = new JFrameMicrographs(fileName);
            frameMicrographs.setLocationRelativeTo(null);
            frameMicrographs.setVisible(true);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            DEBUG.printException(ex);
        }
    }
}
