package xmipp.ij.commons;



public class XmippApplication
{
	
	private static short instances = 0;
    private static short ijwindows = 0;
    private static boolean isscipion = false;
        
	
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
            if(XmippUtil.getXmippImageJ() != null && ijwindows == 0)
            	XmippUtil.getXmippImageJ().setVisible(false);//if ImageJ is not shared hide it
        }
        //System.out.printf("instances %s\n", instances);
		if (instances == 0)
			System.exit(0);
		
	}
        
        public static boolean isScipion()
        {
            return isscipion;
        }
        
        public static void setIsScipion(boolean isscipion)
        {
            XmippApplication.isscipion = isscipion;
        }
        

}
