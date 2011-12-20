/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser;

/**
 *
 * @author Juanjo Vega
 */
public class DEBUG {

    private static boolean DEBUG;

    public static void enableDebug(boolean enable) {
        DEBUG = enable;
    }

    public static void printMessage(String message) {
        if (DEBUG) {
            System.out.println(message);
        }
    }

    public static void printException(Exception ex) {
        if (DEBUG) {
            ex.printStackTrace();
        }
    }
}
