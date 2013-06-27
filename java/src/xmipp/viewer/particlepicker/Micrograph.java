package xmipp.viewer.particlepicker;

import ij.IJ;
import ij.ImagePlus;

import java.awt.Image;
import java.io.File;
import java.util.List;
import java.util.logging.Level;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import xmipp.ij.commons.XmippIJUtil;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.Particle;
import xmipp.utils.XmippMessage;

public abstract class Micrograph {

	private static Icon noimageicon;
	private String file;
	private String name;
	private ImagePlus imp;
	private String pos24file;
	public static final String ext = ".pos";
	public final int width, height;
	private String posfile;
	private String ctf, psd;
	private Icon ctficon;

	public void setPosFileFromXmipp24(String posfile) {
		this.pos24file = posfile;
	}

	public String getPosFileFromXmipp24() {
		return pos24file;
	}

	public Micrograph(String file) {
		this(file, getName(file, 1), null, null);
	}
	
	public Micrograph(String file, String name) {
		this(file, name, null, null);
	}
	
	public Micrograph(String file, String psd, String ctf) {
		this(file, getName(file, 1), psd, ctf);
	}

	public Micrograph(String file, String name, String psd, String ctf) {
		this.file = file;
		this.psd = psd;
		this.ctf = ctf;
		if (!new File(file).exists()) throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("file", file));
		ImageGeneric ig;
		try {
			ig = new ImageGeneric(file);

			width = ig.getXDim();
			height = ig.getYDim();
		} catch (Exception e) {

			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
		this.name = name;
		this.posfile = name + ext;

	}
	

	public ImagePlus getPSDImage()
	{
			if(psd == null || !(new File(psd).exists()))
				return null;
			return XmippIJUtil.getImagePlus(psd);
			
	}
	
	public Icon getCTFIcon()
	{
		String file;
		if(ctficon == null)
		{
			if(psd == null || !(new File(psd).exists()))
				file = (Filename.getXmippPath("resources" + File.separator + "no-image.jpg"));
			else
				file = psd;
			ImagePlus imp = XmippIJUtil.getImagePlus(file);
			Image image = imp.getImage().getScaledInstance(120, 110, Image.SCALE_SMOOTH);
			ctficon = new ImageIcon(image);
			
		}
		return ctficon;
	}
	
	public static Icon getNoImageIcon()
	{
		String file;
		if(noimageicon == null)
		{
			file = (Filename.getXmippPath("resources" + File.separator + "no-image.jpg"));
			ImagePlus imp = XmippIJUtil.getImagePlus(file);
			Image image = imp.getImage().getScaledInstance(120, 110, Image.SCALE_SMOOTH);
			noimageicon = new ImageIcon(image);
			
		}
		return noimageicon;
	}
	
	public String getPSD()
	{
		return psd;
	}
	
	public String getCTF()
	{
		return ctf;
	}

	public boolean fits(int x, int y, int size) {
		if (x < 0 || y < 0) return false;

		int radius = size / 2;
		if (x - radius < 0) return false;
		if (x + radius > width) return false;
		if (y - radius < 0) return false;
		if (y + radius > height) return false;
		return true;
	}

	public static String getName(String file, int level) {
		if(file == null)
			return null;
		// level can start at 1 for file name, 2 is for parent directory name
		String[] tokens = file.split(File.separator);
		if (tokens.length < level)
			throw new IllegalArgumentException(String.format("Name for micrograph is taken from level %s, invalid path ", level, file));
		String name = tokens[tokens.length - level];
		if (level == 1) {
			int pos = name.lastIndexOf('.');
			if (pos != -1) name = name.substring(0, pos);
		}
		return name;
	}

	public String getPosFile() {
		return posfile;
	}

	public ImagePlus getImagePlus() {
		try {

			if (imp == null) {
				imp = XmippImageConverter.loadImage(file);
				if (imp == null) imp = new ImagePlus(file);
			}
			return imp;
		} catch (Exception e) {
			ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public ImagePlus getImagePlus(List<IJCommand> filters) {
		try {
			if (filters.isEmpty()) return getImagePlus();

			if (imp == null) {
				ImageGeneric ig = new ImageGeneric(file);
				ig.read(ImageGeneric.FIRST_IMAGE);

				// Smooth filter should be the first one
				// because it is applied in Xmipp
				for (IJCommand f : filters)
					if (f.getCommand().equals(ParticlePicker.xmippsmoothfilter)) {
						ig.convert2Datatype(ImageGeneric.UChar);
						ImageGeneric igsmooth = new ImageGeneric(ImageGeneric.UChar);
						igsmooth.resize(ig.getXDim(), ig.getYDim());
						ig.smooth(igsmooth);
						break;
					}
				imp = XmippImageConverter.convertToImagePlus(ig);
				ig.destroy();

			}
			return imp;
		} catch (Exception e) {
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void runImageJFilters(List<IJCommand> filters) {
		for (IJCommand f : filters)
			if (!f.getCommand().equals(ParticlePicker.xmippsmoothfilter)) // this filter was applied
				IJ.run(imp, f.getCommand(), f.getOptions());
	}

	public void releaseImage() {
		imp = null;
	}

	public String getFile() {
		return file;
	}

	public String getName() {
		return name;
	}

	public String toString() {
		return name;
	}

	public abstract boolean hasData();


}
