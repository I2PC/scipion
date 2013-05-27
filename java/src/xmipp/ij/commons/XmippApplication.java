package xmipp.ij.commons;

import xmipp.utils.TasksManager;

public class XmippApplication
{
	
	private static short instances = 0;
	
	public static short getInstances()
	{
		return instances;
	}
	
	public static void addInstance()
	{
		instances ++;
		//System.out.printf("instances:%s\n", instances);
	}
	
	public static void removeInstance()
	{
		instances --;
		//System.out.printf("instances:%s\n", instances);
		if (getInstances() == 0)
		{
			//System.out.println("exit");
			TasksManager.getInstance().stop();
			System.exit(0);
		}
	}

}
