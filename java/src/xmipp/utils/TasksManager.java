package xmipp.utils;

import java.util.concurrent.LinkedBlockingQueue;

public class TasksManager
{

	private static TasksManager tm;
	private LinkedBlockingQueue<Task> queue;
	private Thread consumer;
	private boolean stop;

	private TasksManager()
	{
		queue = new LinkedBlockingQueue<Task>();
		consumer = new Thread(new Consumer(queue));
		consumer.start();
	}

	public static TasksManager getInstance()
	{
		if (tm == null)
			tm = new TasksManager();
		return tm;
	}

	public void addTask(Task t)
	{
		try
		{
			queue.put(t);
		}
		catch (InterruptedException e)
		{
			throw new IllegalArgumentException(e);
		}
	}
	
	public void stop()
	{
		stop = true;
	}

	public class Consumer implements Runnable
	{

		private LinkedBlockingQueue<Task> queue;

		public Consumer(LinkedBlockingQueue<Task> queue)
		{
			this.queue = queue;
		}

		@Override
		public void run()
		{
			try
			{
				while (true && !stop)
					queue.take().doTask();
			}
			catch (InterruptedException e)
			{
				throw new IllegalArgumentException(e);
			}

		}

	}

}
