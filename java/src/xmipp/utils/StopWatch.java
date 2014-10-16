/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.utils;

/**
 *
 * @author airen
 */
public class StopWatch {

    /* Private Instance Variables */
    /** Stores the start time when an object of the StopWatch class is initialized. */
    private long startTime;
    private static StopWatch stopwatch;
    /**
     * Custom constructor which initializes the {@link #startTime} parameter.
     */
    private StopWatch() {
        startTime = System.currentTimeMillis();
    }
    
    public static StopWatch getInstance()
    {
        if(stopwatch == null)
            stopwatch = new StopWatch();
        return stopwatch;
    }

    /**
     * Gets the elapsed time (in seconds) since the time the object of StopWatch was initialized.
     * 
     * @return Elapsed time in seconds.
     */
    public double getElapsedTime() {
        long endTime = System.currentTimeMillis();
        return (double) (endTime - startTime) / (1000);
    }
    
    public void printElapsedTime(String comment) {
//        double time = getElapsedTime();
//        if(comment != null)
//            System.out.printf("%.2f %s\n", time, comment);
//        else
//            System.out.println(time);
    }
    
    public void restart()
    {
        startTime = System.currentTimeMillis();
    }
}
