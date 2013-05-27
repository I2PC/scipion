package xmipp.viewer;

import java.util.LinkedList;

import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import ij.IJ;
import xmipp.jni.Filename;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;
import xmipp.viewer.windows.ImagesWindowFactory;

public class Viewer
{
	public static void main(String args[])
	{
		// Schedule a job for the event dispatch thread:
		// creating and showing this application's GUI.
		final String[] myargs = args;
		SwingUtilities.invokeLater(new Runnable()
		{
			public void run()
			{
				
				Param parameters = new Param(myargs);
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

	static void openFiles(Param parameters) throws Exception
	{
		ImagesWindowFactory.openFilesAsDefault(parameters.files, parameters);
	}

	static String[][] getSeparatedFiles(String items[])
	{
		LinkedList<String> images = new LinkedList<String>();
		LinkedList<String> files = new LinkedList<String>();
		LinkedList<String> storeTo;

		for (int i = 0; i < items.length; i++)
		{
			String filename = items[i];

			try
			{
				storeTo = Filename.isSingleImage(filename) ? images : files;
				storeTo.add(filename);
			}
			catch (Exception e)
			{
				IJ.error(e.getMessage());
			}
		}

		String result[][] = new String[2][];

		if (!images.isEmpty())
		{
			result[0] = images.toArray(new String[images.size()]);
		}

		if (!files.isEmpty())
		{
			result[1] = files.toArray(new String[files.size()]);
		}

		return result;
	}
}