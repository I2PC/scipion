/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers;

import browser.ICONS_MANAGER;
import browser.LABELS;
import browser.windows.ImageWindowOperations;
import browser.windows.ImagesWindowFactory;
import browser.windows.StackWindowOperations;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.JButton;

/**
 *
 * @author Juanjo Vega
 */
public class JDialogXmippBrowser extends JDialogXmippFilesList {

    protected JButton jbCaptureWindow;

    public JDialogXmippBrowser(String directory) {
        this(directory, "");
    }

    public JDialogXmippBrowser(String directory, String expression) {
        this(directory, expression, false);
    }

    public JDialogXmippBrowser(String directory, boolean singleSelection) {
        this(directory, "", singleSelection);
    }

    public JDialogXmippBrowser(String directory, String expression, boolean singleSelection) {
        super(directory, 0, expression, singleSelection);   // Port won't be used as method is overriden.

        setModal(false);
        setTitle(LABELS.TITLE_XMIPP_BROWSER);

        // Toolbar buttons.
        jbCaptureWindow = new JButton(LABELS.BUTTON_CAPTURE_WINDOW, ICONS_MANAGER.CAPTURE_WINDOW);
        jbCaptureWindow.setFocusable(false);
        jbCaptureWindow.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbCaptureWindow.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbCaptureWindow.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent ae) {
                captureFrames();
            }
        });

        jToolBar.add(jbCaptureWindow);

        jbOk.setText(LABELS.BUTTON_OPEN_AS_IMAGE);
        jbCancel.setText(LABELS.BUTTON_OPEN_AS_TABLE);
    }

    @Override
    void cancel() {
    }

    @Override
    protected void button1Clicked() {
        panelXmippBrowser.openSelectedFilesAsImages();
    }

    @Override
    protected void button2Clicked() {
        panelXmippBrowser.openSelectedFilesAsTable();
    }

    private void captureFrames() {
        int ids[] = WindowManager.getIDList();
        ArrayList<String> windows = new ArrayList<String>();

        if (ids != null) {
            for (int i = 0; i < ids.length; i++) {
                ImageWindow window = WindowManager.getImage(ids[i]).getWindow();

                if (!(window instanceof ImageWindowOperations) && !(window instanceof StackWindowOperations)) {
                    windows.add(window.getTitle());
                }
            }
        }

        GenericDialog gd = new GenericDialog("Capture frame");
        String titles[] = new String[windows.size()];
        windows.toArray(titles);
        gd.addChoice("Windows:", titles, windows.size() > 0 ? titles[0] : "");
        gd.showDialog();
        if (windows.size() > 0 && !gd.wasCanceled()) {
            String selected = gd.getNextChoice();
            ImagesWindowFactory.captureFrame(WindowManager.getImage(selected));
        }
    }
}
