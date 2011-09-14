package model;

public class Constants {
	

	public static String getIllegalDeleteMsg(String item)
	{
		return String.format("There it must be at least one %s defined", item);

	}
	
	public static String getAlreadyExistsGroupNameMsg(String name)
	{
		return "Group " + name 	+ " already exists";
	}

	public static String getEmptyFieldMsg(String field) {
		return String.format("Must specify a %s", field);
	}

	public static String getAssociatedDataMsg(String field) {
		return field + " has associated data. Can not be removed";
	}

	public static String getNoSuchFieldValueMsg(String field, Object value) {
		return String.format("No such %s %s exists", field, value);
	}
	
	
}
