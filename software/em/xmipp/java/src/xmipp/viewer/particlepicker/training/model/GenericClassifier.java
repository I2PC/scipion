package xmipp.viewer.particlepicker.training.model;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;
import javax.swing.JFrame;
import xmipp.jni.Classifier;
import xmipp.jni.Particle;

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
	
	
	@Override
	public Particle[] autopick(String micrograph, int percent, JFrame frame)
	{
		String command = properties.getProperty("command");
		for(Classifier.Parameter param: params)
		{
			command = command.replace("%(" + param.name + ")", param.value);
		}
		command = command.replace("%(micrograph)", micrograph);
		System.out.println(command);
		try
		{
			//XmippWindowUtil.executeGUICommand(command, frame, "Autopicking ...");
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	@Override
	public int getTrainingParticlesMinimum()
	{
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setParticleSize(int psize)
	{
		// TODO Auto-generated method stub
		
	}


	@Override
	public boolean needsTraining()
	{
		return false;
	}
}
