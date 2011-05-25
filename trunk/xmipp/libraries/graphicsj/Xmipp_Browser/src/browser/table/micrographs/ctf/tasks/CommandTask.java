/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.ctf.tasks;

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
            /*            int total = 5000;
            int delay = 1000;
            while (total > 0) {
            System.out.println(total + " seconds to finish: "
            + command.substring(0, 20) + (row >= 0 ? "/ row=" + row : ""));
            Thread.sleep(delay);
            total -= delay;
            }*/
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
        System.out.println(">> " + command);
        Process p = Runtime.getRuntime().exec(command);

        BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

        // read the output from the command
        String s;
        while ((s = stdInput.readLine()) != null) {
            System.out.println(s);
        }
    }
}
