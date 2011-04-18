package browser.table.micrographs.ctf.tasks;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class TasksEngine implements iCommandsListener {

    private int tasks = 0;
    private iCommandsListener commandsListener;

    public TasksEngine(iCommandsListener commandsListener) {
        this.commandsListener = commandsListener;
    }

    public synchronized void add(CommandTask task) {
        tasks++;
        task.start();
    }

    public synchronized boolean isDone() {
        return tasks == 0;
    }

    public synchronized void done() {
        tasks--;

        if (isDone()) {
            commandsListener.done();
        }
    }
}
