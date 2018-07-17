package xmipp.viewer.particlepicker.training.model;

public enum Mode {
	Manual, // Only user-selected particles
	Supervised, // User trains the picker, some particles are manual and others automatic
    Available, // Micrograph is available to autopick in Xmipp
    Review, // User checks results from autopicking and can edit
    ReadOnly,  // User can only see the coordinates
    Extract,
    Automatic; // Used in GenericPicker, user can setup parameters, but not modify coordinates manually
	
	
	public static Mode getMode(String s)
	{
		if(s.equalsIgnoreCase(Mode.Manual.toString()))
			return Mode.Manual;
		else if(s.equalsIgnoreCase(Mode.Supervised.toString()))
			return Mode.Supervised;
		else if(s.equalsIgnoreCase(Mode.Review.toString()))
			return Mode.Review;
		else if(s.equalsIgnoreCase(Mode.Available.toString()))
			return Mode.Available;
		else if(s.equalsIgnoreCase(Mode.ReadOnly.toString()))
			return Mode.ReadOnly;
		else if(s.equalsIgnoreCase(Mode.Extract.toString()))
			return Mode.Extract;
		else if (s.equalsIgnoreCase(Mode.Automatic.toString()))
		    return Mode.Automatic;

		return null;
	}

}
