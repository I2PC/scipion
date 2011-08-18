package ppicker;

public class Constants {
	
	public static String empty_group_name_msg = "Must specify a Group Name";
	public static String delete_last_group_msg = "There it must be at least one group defined.";

	public static String getAlreadyExistsGroupNameMsg(String name)
	{
		return "Group " + name 	+ " already exists";
	}
}
