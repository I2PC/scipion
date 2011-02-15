/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tests;

import java.util.TimerTask;

/**
 *
 * @author Juanjo Vega
 */
public class Abort extends TimerTask {

    public Abort() {
    }

    @Override
    public void run() {
        System.err.println(" *** Aborting...");
        System.exit(1);
    }
}
