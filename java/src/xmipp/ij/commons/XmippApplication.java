package xmipp.ij.commons;

import ij.WindowManager;


public class XmippApplication
{
	
	private static short instances = 0;
        private static short ijwindows = 0;
	
	public static short getInstances()
	{
		return instances;
	}
	
	public static void addInstance(boolean isijwindow)
	{
		instances ++;
                if(isijwindow)
                    ijwindows ++;
		//System.out.printf("instances:%s\n", instances);
	}
        
        
	
	public static void removeInstance(boolean isijwindow)
	{
		instances --;
                if(isijwindow)
                {
                    ijwindows --;
                    if(XmippIJUtil.getXmippImageJ() != null && ijwindows == 0)
			XmippIJUtil.getXmippImageJ().setVisible(false);//if ImageJ is not shared hide it
                }
                //System.out.printf("instances %s\n", instances);
		if (instances == 0)
			System.exit(0);
		
	}

}
