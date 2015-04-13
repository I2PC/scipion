package xmipp.viewer.particlepicker.training.model;
public enum Mode {
	Manual, Supervised, Available, Review, ReadOnly, Extract;
	
	
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
		return null;
	}

}
