/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.micrographs.ctf.tasks;

import xmipp.utils.DEBUG;
import ij.IJ;
import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 *
 * @author Juanjo Vega
 */
public class CommandTask implements Runnable {

    protected String command;
    protected int row = -1;
    protected iTaskCompletionListener commandsListener;

    public CommandTask(String command, int row, iTaskCompletionListener commandsListener) {
        this(command, commandsListener);
        this.row = row;
    }

    public CommandTask(String command, iTaskCompletionListener commandsListener) {
        this(command);
        this.commandsListener = commandsListener;
    }

    public CommandTask(String command) {
        this.command = command;
    }

    public void start() {
        Thread t = new Thread(this);
        t.start();
    }

    public void run() {
        try {
            runTask();
        } catch (Exception ex) {
            IJ.error("Error running command: " + ex.getMessage());
        } finally {
            if (commandsListener != null) {
                commandsListener.done(this);
            }
        }
    }

    private void runTask() throws Exception {
        DEBUG.printMessage(">> " + command);
        Process p = Runtime.getRuntime().exec(command);

        BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

        // read the output from the command
        String s;
        while ((s = stdInput.readLine()) != null) {
            DEBUG.printMessage(s);
        }
    }
}
