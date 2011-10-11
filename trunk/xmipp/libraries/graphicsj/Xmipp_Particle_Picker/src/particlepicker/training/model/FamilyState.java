package particlepicker.training.model;
public enum FamilyState {
	Manual, Supervised, Available, Review;
	
	
	public static FamilyState getFamilyState(String s)
	{
		if(s.equalsIgnoreCase(FamilyState.Manual.toString()))
			return FamilyState.Manual;
		if(s.equalsIgnoreCase(FamilyState.Supervised.toString()))
			return FamilyState.Supervised;
		if(s.equalsIgnoreCase(FamilyState.Review.toString()))
			return FamilyState.Review;
		if(s.equalsIgnoreCase(FamilyState.Available.toString()))
			return FamilyState.Available;
		return null;
	}

}
