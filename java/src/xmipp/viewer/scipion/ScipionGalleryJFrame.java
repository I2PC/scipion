/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import xmipp.jni.MetaData;
import xmipp.utils.Param;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryJFrame extends GalleryJFrame {

    private final String cmdname;
    private final String cmdscript;

    public ScipionGalleryJFrame(String filename, MetaData md, Param parameters) {
        super(filename, md, parameters);
        cmdname = parameters.cmdname;
        cmdscript = parameters.cmdscript;
        initComponents();
    }

    private void initComponents() {
        if (cmdname != null) {
            JButton cmdbutton = XmippWindowUtil.getTextButton(cmdname, new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    runScipionScriptOnSelection();
                }
            });

            commandspn.add(cmdbutton);
        }
    }

    public void runScipionScriptOnSelection() {
        try {

            String selectionmd = "selection.xmd";
            saveSelection(selectionmd, true);
            String pwhome = System.getenv("PW_HOME");
            String command = String.format("%s/pw.bashrc\n %s/apps/%s %s", pwhome, pwhome, cmdscript, selectionmd);
            System.out.println(command);
            Process exec = Runtime.getRuntime().exec(command);
            exec.getou
        } catch (Exception ex) {

            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
