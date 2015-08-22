package xmipp.ij.commons;

public class IJCommand
{

	private String command;
	private String options;
	
	public IJCommand(String command, String options)
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
	
    public String getMdCommand()
    {
        return command.replace('_', '&').replace(' ', '_');
    }
    
    public String getMdOptions()
    {
        return (options == null || options.equals("")) ? "NULL" :options.replace('_', '&').replace(' ', '_');
    }

    public static String getString(String mdstring)
    {
        return mdstring.replace('_', ' ').replace('&', '_');
    }
	
    
    public String toString()
    {
    	return command + " " + options;
    }
	
}
