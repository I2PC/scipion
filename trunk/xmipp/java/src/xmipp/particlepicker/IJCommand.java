package xmipp.particlepicker;

public class Filter
{

	private String command;
	private String options;
	
	public Filter(String command, String options)
	{
		this.command = command;
		this.options = options;
	}
	
	public String getCommand()
	{
		return command;
	}
	public void setCommand(String command)
	{
		this.command = command;
	}
	public String getOptions()
	{
		return options;
	}
	public void setOptions(String options)
	{
		this.options = options;
	}
	
	
	
	
	
}
