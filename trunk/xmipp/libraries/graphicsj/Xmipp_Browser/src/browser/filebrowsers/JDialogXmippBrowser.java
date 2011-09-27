/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers;

import browser.COMMAND_PARAMETERS;
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
import javax.swing.JToolBar;

/**
 *
 * @author Juanjo Vega
 */
public class JDialogXmippBrowser extends JDialogXmippFilesList {

    JButton jbParent, jbRefresh, jbCaptureWindow;
    JToolBar jToolBar;

    public JDialogXmippBrowser(String directory) {
        this(directory, "");
    }

    public JDialogXmippBrowser(String directory, String expression) {
        this(directory, false, COMMAND_PARAMETERS.SELECTION_TYPE_ANY, expression);
    }

    public JDialogXmippBrowser(String directory, boolean singleSelection) {
        this(directory, singleSelection, COMMAND_PARAMETERS.SELECTION_TYPE_ANY, "");
    }

    public JDialogXmippBrowser(String directory, boolean singleSelection, String expression) {
        this(directory, singleSelection, COMMAND_PARAMETERS.SELECTION_TYPE_ANY, expression);
    }

    public JDialogXmippBrowser(String directory, boolean singleSelection, String selType, String expression) {
        super(directory, 0, singleSelection, selType, expression);   // Port won't be used as method is overriden.

//        setModal(false);
//        setAlwaysOnTop(false);
        toFront();
        setTitle(LABELS.TITLE_XMIPP_BROWSER);

        setToolbar();

        jbOk.setText(LABELS.BUTTON_OPEN_AS_IMAGE);
        jbCancel.setText(LABELS.BUTTON_OPEN_AS_GALLERY);
    }

    void setToolbar() {
        jToolBar = new javax.swing.JToolBar();
        jbParent = new javax.swing.JButton();
        jbRefresh = new javax.swing.JButton();

        jToolBar.setFloatable(false);
        jToolBar.setRollover(true);

        jbParent.setIcon(ICONS_MANAGER.PARENT_DIRECTORY);
        jbParent.setText(LABELS.BUTTON_PARENT_DIRECTORY);
        jbParent.setFocusable(false);
        jbParent.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbParent.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbParent.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                goParent();
            }
        });

        jbRefresh.setIcon(ICONS_MANAGER.REFRESH_DIRECTORY);
        jbRefresh.setText(LABELS.BUTTON_REFRESH_DIRECTORY);
        jbRefresh.setFocusable(false);
        jbRefresh.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbRefresh.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbRefresh.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                refresh();
            }
        });

        jbCaptureWindow = new JButton(LABELS.BUTTON_CAPTURE_WINDOW, ICONS_MANAGER.CAPTURE_WINDOW);
        jbCaptureWindow.setFocusable(false);
        jbCaptureWindow.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbCaptureWindow.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbCaptureWindow.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent ae) {
                captureFrames();
            }
        });

        jToolBar.add(jbParent);
        jToolBar.add(jbRefresh);
        jToolBar.add(jbCaptureWindow);

        getContentPane().add(jToolBar, java.awt.BorderLayout.NORTH);
    }

    void refresh() {
        panelXmippBrowser.refreshCurrentDirectory();
    }

    protected void goParent() {
        panelXmippBrowser.goParent();
    }

    @Override
    protected void cancel() {
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
