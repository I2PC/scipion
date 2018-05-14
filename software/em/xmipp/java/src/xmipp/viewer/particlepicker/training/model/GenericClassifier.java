package xmipp.viewer.particlepicker.training.model;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;

import javax.swing.JFrame;

import xmipp.jni.Classifier;
import xmipp.jni.MetaData;
import xmipp.jni.Particle;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.DEBUG;


public class GenericClassifier extends Classifier
{

	protected String classifierProperties;
	protected Properties properties;
	

	public GenericClassifier(String classifierProperties)
	{
		try
		{
			this.classifierProperties = classifierProperties;
			// create and load properties
			properties = new Properties();
			FileInputStream in = new FileInputStream(classifierProperties);
			properties.load(in);
			
			in.close();
			properties.setProperty("applyChanges", String.valueOf(applyChanges));
			String[] paramNames = properties.getProperty("parameters").split(",");
			Parameter param;
			String label, help, value;
			params = new ArrayList<Classifier.Parameter>();
			for(String paramName : paramNames)
			{
				label = properties.getProperty(paramName + ".label");
				help = properties.getProperty(paramName + ".help");
				value = properties.getProperty(paramName + ".value");
				param = new Parameter(paramName, label, help, value);
				params.add(param);
			}
			if (hasInitialCoordinates())
				convert();
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	/** Execute a command and print it if in DEBUG mode */
	private void execute(String name, String cmd, String runDir) throws Exception
	{
		if(cmd == null)
		{
			if (DEBUG.hasScipionDebug())
				System.out.println(name + " command: NULL\n");
			return;
		}

		if (DEBUG.hasScipionDebug())
			System.out.println(name + " command: \n" + cmd + "\n");

		String output = (runDir == null) ? XmippWindowUtil.executeCommand(cmd, true) :
				                           XmippWindowUtil.executeCommand(cmd, true, runDir);
		if (DEBUG.hasScipionDebug())
			System.out.println(name + " output: \n" + output + "\n");
	}

	private void convert() throws Exception {
		String convertCommand = properties.getProperty("convertCommand");
		execute("Convert", convertCommand, null);
	}

	public void autopick(SupervisedPickerMicrograph micrograph)
	{
		String preprocessCommand = properties.getProperty("preprocessCommand");
		String autopickCommand = properties.getProperty("autopickCommand");
		String runDir = properties.getProperty("runDir");

		for(Classifier.Parameter param: params)
		{
			if(preprocessCommand != null)
				preprocessCommand = preprocessCommand.replace("%(" + param.name + ")", param.value);
			autopickCommand = autopickCommand.replace("%(" + param.name + ")", param.value);
		}

		String micPath = micrograph.getFile();

		if (runDir != null)
			micPath = new File(runDir).toURI().relativize(new File(micPath).toURI()).getPath();

		String name = micrograph.getName();
		autopickCommand = autopickCommand.replace("%(micrograph)", micPath);
		autopickCommand = autopickCommand.replace("%(micrographName)", name);

		try
		{
			execute("Preprocess", preprocessCommand, runDir);
			execute("Autopick", autopickCommand, runDir);
			convert();
			writeProperties();
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	public int getTrainingParticlesMinimum()
	{
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setSize(int psize)
	{
		// TODO Auto-generated method stub
	}

	@Override
	public boolean needsTraining()
	{
		return false;
	}
	
	public void writeProperties()
	{
		try
		{
			for(Classifier.Parameter param: params)
				properties.setProperty(param.name + ".value", param.value);
			properties.setProperty("applyChanges", String.valueOf(applyChanges));
			File file = new File(classifierProperties);
			FileOutputStream fileOut;
			fileOut = new FileOutputStream(file);
			properties.store(fileOut, "Favorite Things");
			fileOut.close();
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
			if(!e.getMessage().isEmpty())
				XmippDialog.showError(null,e.getMessage());
		}
	}
	
	public void setApplyChanges(boolean applyChanges) {
		super.setApplyChanges(applyChanges);
		writeProperties();
	}

	public boolean hasInitialCoordinates() {
		return properties.getProperty("hasInitialCoordinates", "").equals("true");
	}

	public boolean doPickAll() {
		return properties.getProperty("doPickAll", "").equals("true");
	}

}
