package xmipp.jni;

import java.util.List;

import javax.swing.JFrame;

public abstract class Classifier
{
	protected List<Classifier.Parameter> params;
	
	public abstract int getTrainingParticlesMinimum();
	
	public abstract void setParticleSize(int psize);
	
	public abstract Particle[] autopick(String micrograph, int percent, JFrame frame);
	
	public class Parameter
	{
		public String name;
		public String label;
		public String help;
		public String value;
		
		public Parameter(String name, String label, String help, String value)
		{
			this.name = name;
			this.label = label;
			this.help = help;
			this.value = value;
		}
		
	}
	
	public List<Classifier.Parameter> getParams()
	{
		return params;
	}
	
	public abstract boolean needsTraining();


}
