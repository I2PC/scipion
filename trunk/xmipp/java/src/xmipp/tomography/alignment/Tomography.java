package xmipp.tomography.alignment;

public class Tomography {
	
	private int tiltangle;
	private String tomofile;
	
	public Tomography(String tomofile, int tiltangle)
	{
		this.tomofile = tomofile;
		this.tiltangle = tiltangle;
	}

}
