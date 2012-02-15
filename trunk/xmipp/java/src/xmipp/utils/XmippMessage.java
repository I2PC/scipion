package xmipp.utils;


/** This class will contains methods to return formatted text messages */
public class XmippMessage {
	

	/** Particle Picker messages */
	public static String getIllegalDeleteMsg(String item)
	{
		return String.format("There it must be at least one %s defined", item);
	}
	
	public static String getAlreadyExistsGroupNameMsg(String name)
	{
		return String.format("Group %s already exists", name);
	}

	public static String getEmptyFieldMsg(String field) {
		return String.format("Must specify %s", field);
	}

	public static String getAssociatedDataMsg(String field) {
		return String.format("%s has associated data. Can not be removed");
	}

	public static String getNoSuchFieldValueMsg(String field, Object value) {
		return String.format("No such %s %s exists", field, value);
	}

	public static String getOutOfBoundsMsg(Object o)
	{
		return String.format("%s out of bounds", o);
	}
	
	
}
