/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.utils;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.Timer;

/**
 *
 * @author Juanjo Vega
 */
public class TaskTimer implements ActionListener {

    Runnable runnable;
    Timer timer;

    public TaskTimer(Runnable runnable) {
        this.runnable = runnable;
        timer = new Timer(500, this);
        timer.setRepeats(false);    // Repeats only once.
    }

    public void start() {
        if (timer.isRunning()) {
            timer.restart();
        } else {
            timer.start();
        }
    }

    public boolean isRunning() {
        return timer.isRunning();
    }

    public void actionPerformed(ActionEvent e) {
        runnable.run();
    }
}
