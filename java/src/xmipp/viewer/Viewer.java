package xmipp.viewer;

import javax.swing.SwingUtilities;
import xmipp.utils.DEBUG;
import xmipp.utils.Params;
import xmipp.utils.StopWatch;

import xmipp.viewer.windows.ImagesWindowFactory;

public class Viewer
{
	public static void main(String args[])
	{
                StopWatch stopWatch = StopWatch.getInstance();
                stopWatch.printElapsedTime("starting app");
		// Schedule a job for the event dispatch thread:
		// creating and showing this application's GUI.
		final String[] myargs = args;
		SwingUtilities.invokeLater(new Runnable()
		{
			public void run()
			{
				Params parameters = new Params(myargs);
				try
				{
					if (parameters.debug)
						DEBUG.enableDebug(true);

					if (parameters.files != null)
					{
						// DEBUG.enableDebug(true);
						openFiles(parameters);
					}
				}
				catch (Exception e)
				{
					e.printStackTrace();
				}

			}
		});
	}
        
  

	static void openFiles(Params parameters) throws Exception
	{
		ImagesWindowFactory.openFilesAsDefault(parameters.files, parameters);
	}

	
}