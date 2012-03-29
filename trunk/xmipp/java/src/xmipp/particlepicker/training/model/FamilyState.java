package xmipp.particlepicker.training.model;
public enum FamilyState {
	Manual, Supervised, Available, Review, ReadOnly;
	
	
	public static FamilyState getFamilyState(String s)
	{
		if(s.equalsIgnoreCase(FamilyState.Manual.toString()))
			return FamilyState.Manual;
		else if(s.equalsIgnoreCase(FamilyState.Supervised.toString()))
			return FamilyState.Supervised;
		else if(s.equalsIgnoreCase(FamilyState.Review.toString()))
			return FamilyState.Review;
		else if(s.equalsIgnoreCase(FamilyState.Available.toString()))
			return FamilyState.Available;
		else if(s.equalsIgnoreCase(FamilyState.ReadOnly.toString()))
			return FamilyState.ReadOnly;
		return null;
	}

}
