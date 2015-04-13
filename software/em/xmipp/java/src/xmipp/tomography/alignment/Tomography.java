package xmipp.tomography.alignment;


import ij.ImagePlus;

import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.labun.surf.Descriptor;
import com.labun.surf.Detector;
import com.labun.surf.IntegralImage;
import com.labun.surf.InterestPoint;
import com.labun.surf.Matcher;
import com.labun.surf.Params;

import xmipp.ij.commons.XmippUtil;

public class Tomography {
	
	private int tiltangle;
	private String tomofile;
	private Tomography previous;
	private Tomography next;
	private Matrix atmatrix;
	private ImagePlus imp;
	private IntegralImage intimg;
	
	public Tomography(String tomofile, int tiltangle)
	{
		this.tomofile = tomofile;
		this.tiltangle = tiltangle;
	}
	
	public Tomography(String tomofile, int tiltangle, Tomography previous) {
		this(tomofile, tiltangle);
		setPrevious(previous);
	}

	
	public void setPrevious(Tomography previous)
	{
		this.previous = previous;
	}
	
	public void setNext(Tomography next)
	{
		this.next = next;
	}

	
	public Tomography getPrevious() 
	{
		return previous;
	}
	
	public Tomography getNext()
	{
		return next;
	}
	
	public void computeAffineTransform()
	{
		
		Map<InterestPoint, InterestPoint> cpm = getCommonPointsMapping();
		if (cpm == null || cpm.isEmpty())
			return;//no affine transform
		System.out.println("Computing affine transform on  " + tomofile + " and " + getPrevious().getTomofile());
		
		System.out.println(cpm.size());
	}
	
	
	private String getTomofile() {
		return tomofile;
	}

	public Map<InterestPoint, InterestPoint> getCommonPointsMapping()
	{
		if(previous == null)
			return null;
		
		Params pparams = new Params(4, 4, 0.0001f, 2, false, false, false, 1, false);//octaves, layers, threshold, initialStep, upright, ...
		Params params = new Params(pparams);
		
		
		List<InterestPoint> pipoints = detectAndDescribeInterestPoints(previous.getIntegralImage(), pparams);
		List<InterestPoint> ipoints = detectAndDescribeInterestPoints(getIntegralImage(), params);
		
		
		Map<InterestPoint, InterestPoint> matchedpoints = Matcher.findMathes(pipoints, ipoints);		
		Map<InterestPoint, InterestPoint> revmatchedpoints = Matcher.findMathes(ipoints, pipoints);
		// take only those points that matched in the reverse comparison too
		Map<InterestPoint, InterestPoint> commonpoints = new HashMap<InterestPoint, InterestPoint>();
		for (InterestPoint ppoint : matchedpoints.keySet()) {
			InterestPoint point = matchedpoints.get(ppoint);
			if (ppoint == revmatchedpoints.get(point))
				commonpoints.put(ppoint, point);
		}
		
		return commonpoints;
	}
		
	public static List<InterestPoint> detectAndDescribeInterestPoints(IntegralImage intImg, Params p) {
		
		List<InterestPoint> ipts = Detector.fastHessian(intImg, p);
		
		
		// Describe interest points with SURF-descriptor
		if (!p.isUpright())
			for (InterestPoint ipt: ipts)
				Descriptor.computeAndSetOrientation(ipt, intImg);
		for (InterestPoint ipt: ipts)
			Descriptor.computeAndSetDescriptor(ipt, intImg, p);

		return ipts;
	}
	
	public Matrix getAffineTransform()
	{
		return atmatrix;
	}
	
	public Rectangle getOverlappedArea()
	{
		return null;
	}

	public ImagePlus getImagePlus()
	{
		if(imp == null)
			imp = XmippUtil.getImagePlus(tomofile);
		return imp;
	}
	
	public IntegralImage getIntegralImage()
	{
		if(intimg == null)
			intimg = new IntegralImage(getImagePlus().getProcessor(), true);
		return intimg;
	}
	
	
	

}
