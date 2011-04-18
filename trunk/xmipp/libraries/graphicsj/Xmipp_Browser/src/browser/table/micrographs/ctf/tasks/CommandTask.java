/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.micrographs.ctf.tasks;

import ij.IJ;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 *
 * @author Juanjo Vega
 */
public class CommandTask implements Runnable {

    protected String command;
    protected iCommandsListener commandsListener;

    public CommandTask(String command, iCommandsListener commandsListener) {
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
            System.out.println(">> " + command);
            Process p = Runtime.getRuntime().exec(command);

            BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

            // read the output from the command
            String s;
            while ((s = stdInput.readLine()) != null) {
                System.out.println(s);
            }
        } catch (IOException ex) {
            IJ.error("Error running command: " + ex.getMessage());
            throw new RuntimeException(ex);
        } finally {
            if (commandsListener != null) {
                commandsListener.done();
            }
        }
    }
}
