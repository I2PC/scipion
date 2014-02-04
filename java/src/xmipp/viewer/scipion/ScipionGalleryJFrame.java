/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
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

    private final String set;
    private final String script;
    private final String imagesid;

    public ScipionGalleryJFrame(String filename, MetaData md, ScipionParams parameters) {
        super(filename, md, parameters);
        set = parameters.set;
        script = parameters.script;
        imagesid = parameters.imagesid;
        initComponents();
    }

    private void initComponents() {
        if (set != null) {
            JButton cmdbutton = XmippWindowUtil.getTextButton("Create New " + set, new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    try {
                        String selectionmd = "selection.xmd";
                        selectionmd = new File(selectionmd).getAbsolutePath();
                        saveSelection(selectionmd, true);
                        String command = String.format("%s %s %s %s", script, selectionmd, set.replaceAll("\\s",""), imagesid);
                        executeCommand(command);
                    } catch (Exception ex) {
                        Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            });

            commandspn.add(cmdbutton);
        }
    }



    private void executeCommand(String command) throws Exception {

        System.out.println(command);
        StringBuffer output = new StringBuffer();

        Process p;

        p = Runtime.getRuntime().exec(command);
        p.waitFor();
        BufferedReader reader
                = new BufferedReader(new InputStreamReader(p.getInputStream()));

        output.append("Output\n");
        String line = "";
        while ((line = reader.readLine()) != null) {
            output.append(line + "\n");
        }
        reader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
        output.append("Error\n");
        while ((line = reader.readLine()) != null) {
            output.append(line + "\n");
        }
        System.out.println(output.toString());

    }
}
