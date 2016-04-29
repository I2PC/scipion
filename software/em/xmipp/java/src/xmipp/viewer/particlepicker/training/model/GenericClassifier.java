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
		}
		catch (IOException e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	
	public void autopick(SupervisedPickerMicrograph micrograph)
	{
		String preprocessCommand = properties.getProperty("preprocessCommand");
		String autopickCommand = properties.getProperty("autopickCommand");
		String convertCommand = properties.getProperty("convertCommand");
		String runDir = properties.getProperty("runDir");
		for(Classifier.Parameter param: params)
		{
			if(preprocessCommand != null)
				preprocessCommand = preprocessCommand.replace("%(" + param.name + ")", param.value);
			autopickCommand = autopickCommand.replace("%(" + param.name + ")", param.value);
		}
		String micPath = micrograph.getFile();
		if(runDir != null)
			micPath = new File(runDir).toURI().relativize(new File(micPath).toURI()).getPath();
		String name = micrograph.getName();
		autopickCommand = autopickCommand.replace("%(micrograph)", micPath);
		autopickCommand = autopickCommand.replace("%(micrographName)", name);
		String output;
		try
		{
			if(runDir != null)
			{
				if(preprocessCommand != null)
				{
					output = XmippWindowUtil.executeCommand(preprocessCommand, true, runDir);
					//System.out.println("preprocess output \n" + output);
				}
				output = XmippWindowUtil.executeCommand(autopickCommand, true, runDir);
				if (DEBUG.hasScipionDebug()) {
				    System.out.println("Autopick command: \n" + autopickCommand + "\n");
				    System.out.println("Autopick output: \n" + output);
				}
			}
			else
			{
				if(preprocessCommand != null)
				{
					output = XmippWindowUtil.executeCommand(preprocessCommand, true);
					//System.out.println("preprocess output \n" + output);
				}
				output = XmippWindowUtil.executeCommand(autopickCommand, true);
				if (DEBUG.hasScipionDebug()) {
				    System.out.println("Autopick command: \n" + autopickCommand + "\n");
				    System.out.println("Autopick output: \n" + output);
				}
			}
			
			output = XmippWindowUtil.executeCommand(convertCommand, true);
			//System.out.println("convert output \n" + output);
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
}
