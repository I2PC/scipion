package xmipp.utils;

import java.util.ArrayList;
import java.util.List;

public class TasksManager
{

	private List<Runnable> tasks;
	private static TasksManager tm;

	public static TasksManager getInstance()
	{
		if (tm == null)
			tm = new TasksManager();
		return tm;
	}

	private TasksManager()
	{
		tasks = new ArrayList<Runnable>();
	}

	public void addTask(Runnable r)
	{
		tasks.add(r);
		new Thread(r).start();
	}

	

	public void removeTask(Runnable r)
	{
		tasks.remove(r);
	}

	
}
