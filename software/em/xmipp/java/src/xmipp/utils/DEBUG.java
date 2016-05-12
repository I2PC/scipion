/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.utils;
/**
 *
 * @author Juanjo Vega
 */
public class DEBUG {

    private static boolean DEBUG = false;

    public static void enableDebug(boolean enable) {
        DEBUG = enable;
    }

    public static void printMessage(String message) {
        if (DEBUG) {
            System.out.println(message);
        }
    }
    
    public static void printFormat(String format, Object...objects){
    	if (DEBUG){
    		System.out.format(format+"\n", objects);
    	}
    }

    public static void printException(Exception ex) {
        if (DEBUG) {
            ex.printStackTrace();
        }
    }
    
    public static void printStackTrace() {
    	if (DEBUG){
    		Exception ex = new Exception();
    		ex.printStackTrace();
    	}
    }

    public static boolean hasScipionDebug() {
        String debug = System.getenv("SCIPION_DEBUG");
        return (debug != null && debug.equals("1"));
    }

    //Timing
    public static long millisecs;
    
    public static long tic(){
    	millisecs = System.currentTimeMillis();
    	return millisecs;
    }
    
    public static long toc(){
    	long current = System.currentTimeMillis();
    	long result = current - millisecs;
    	millisecs = current;
    	System.out.format("took: %d milliseconds %n", result);
    	return result;
    }

}
