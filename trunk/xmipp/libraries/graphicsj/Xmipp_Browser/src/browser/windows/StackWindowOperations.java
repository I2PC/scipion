package browser.windows;

import browser.DEBUG;
import browser.windows.menubar.XmippMenuBar;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageLayout;
import ij.gui.StackWindow;
import ij.io.FileInfo;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Timer;
import java.util.TimerTask;
import javax.swing.JPanel;
import xmippij.XmippImageConverter;

/**
 *
 * @author Juanjo Vega
 */
public class StackWindowOperations extends StackWindow implements iPollImageWindow {

    protected boolean poll = false;
    protected Timer timer = null;
    protected final static int PERIOD = 5000;  // repeat every 5 seconds
    protected XmippMenuBar menuBar;
    //private String path;// = imp.getOriginalFileInfo().directory + imp.getOriginalFileInfo().fileName;
    private File f;// = new File(path);
    private long last;// = f.lastModified();

    public StackWindowOperations(ImagePlus imp) {
        this(imp, false);
    }

    public StackWindowOperations(ImagePlus imp, boolean poll) {
        super(imp);

        JPanel previousContent = new JPanel();
        previousContent.setLayout(new ImageLayout(ic));

        Component components[] = getComponents();
        for (int i = 0; i < components.length; i++) {
            previousContent.add(components[i]);
        }

        removeAll();
        setLayout(new BorderLayout());

        add(previousContent, BorderLayout.CENTER);

        menuBar = new XmippMenuBar(this, imp);
        setMenuBar(menuBar);

        fixSize();

        setInitialPoll(poll);
    }

    @Override
    public void windowOpened(WindowEvent e) {
        pack();
    }

    public boolean isPoll() {
        return poll;
    }

    protected void setInitialPoll(boolean poll) {
        // Sets if image can poll (reloaded from disk) or not.
        if (menuBar.allowsPolling()) {
            FileInfo ofi = imp.getOriginalFileInfo();
            f = new File(ofi.directory + ofi.fileName);
            last = f.lastModified();

            menuBar.setPollStatus(poll);

            setPoll(poll);  // If polling, starts timer to reload image form disk every period.
        }
    }

    public synchronized void setPoll(boolean poll) {
        System.out.println("setPoll!");
        //     if (this.poll != poll) {    // If state changes.
        this.poll = poll;

        if (poll) {
            startTimer();
        } else {
            if (timer != null) {
                timer.cancel();
                timer = null;
            }
        }
        //   }
    }

    protected void revert() {
        // Reverts only when image has changed since last time.
        if (last < f.lastModified()) {
            DEBUG.printMessage(" *** Reverting from disk...");
            IJ.showStatus("Reloading " + imp.getTitle());

            try {
                XmippImageConverter.revert(imp, f.getAbsolutePath());
            } catch (Exception ex) {
                ex.printStackTrace();
            }

            last = f.lastModified();

            IJ.showStatus("");
        }
    }

    private void startTimer() {
        // Avoid multiple timers running simultaneously.
        if (timer == null) {
            timer = new Timer();

            timer.scheduleAtFixedRate(new TimerTask() {

                public void run() {
                    revert();
                }
            }, 0, PERIOD);
        }
    }

    @Override
    public void windowClosing(WindowEvent e) {
        setPoll(false);

        super.windowClosing(e);
    }

    private void fixSize() {
        pack();

        Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
        Point loc = getLocation();
        Dimension size = getSize();
        if (loc.y + size.height > screen.height) {
            getCanvas().zoomOut(0, 0);
        }
    }
}
