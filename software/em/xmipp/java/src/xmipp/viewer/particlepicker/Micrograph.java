package xmipp.viewer.particlepicker;

import ij.IJ;
import ij.ImagePlus;

import java.awt.Image;
import java.io.File;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import xmipp.ij.commons.IJCommand;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippUtil;
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
	private int width, height;
	private String posfile;
	private CtfInfo ctfInfo;
	private Icon ctficon;

	public void setPosFileFromXmipp24(String posfile) {
		this.pos24file = posfile;
	}

	public String getPosFileFromXmipp24() {
		return pos24file;
	}

	public Micrograph(String file) {
		this(file, getName(file, 1), null);
	}

	public Micrograph(String file, String name) {
		this(file, name, null);
	}

	public Micrograph(String file, CtfInfo ctfInfo) {
		this(file, getName(file, 1), ctfInfo);
	}

	public Micrograph(String file, String name, CtfInfo ctfInfo) {
		this.file = file;
		this.ctfInfo = ctfInfo;
		String path = file;
		if(Filename.hasPrefix(file))
		    path = Filename.removePrefix(file);
		if (!new File(path).exists()) 
        {
             System.out.printf("path %s selfile %s\n", path, ParticlePicker.getPicker().selfile);
             file = Filename.findImagePath(path, ParticlePicker.getPicker().selfile, true);
             if(file == null)
                throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("micrograph", path));
        }
		
		this.name = name;
		this.posfile = name + ext;
                //loadDimensions(); // ensure width and height get updated

	}

	public String getPSD() {
		return ctfInfo == null ? null : ctfInfo.psd;
	}

	public ImagePlus getPSDImage() {
		String psd = getPSD();

		if (psd == null || !(new File(psd).exists()))
			return null;
		return XmippUtil.getImagePlus(psd);

	}

	public Icon getCTFIcon() {
		String file;
		if (ctficon == null) {
			String psd = getPSD();
			if (psd == null || !(new File(psd).exists()))
				file = (Filename.getXmippPath("resources" + File.separator
						+ "no-image.jpg"));
			else
				file = psd;
			ImagePlus imp = XmippUtil.getImagePlus(file);
			Image image = imp.getImage().getScaledInstance(120, 110,
					Image.SCALE_SMOOTH);
			ctficon = new ImageIcon(image);

		}
		return ctficon;
	}

	public static Icon getNoImageIcon() {
		String file;
		if (noimageicon == null) {
			file = (Filename.getXmippPath("resources" + File.separator
					+ "no-image.jpg"));
			ImagePlus imp = XmippUtil.getImagePlus(file);
			Image image = imp.getImage().getScaledInstance(120, 110,
					Image.SCALE_SMOOTH);
			noimageicon = new ImageIcon(image);

		}
		return noimageicon;
	}

	public String getCTF() {
		return ctfInfo == null ? null : ctfInfo.ctf;
	}

    public CtfInfo getCtfInfo() {
        return ctfInfo;
    }

    public boolean fits(int x, int y, int size) {
		if (x < 0 || y < 0)
			return false;

		int radius = size / 2;
		if (x - radius < 0 || x + radius > getWidth() || y - radius < 0
				|| y + radius > getHeigth())
			return false;
		return true;
	}

	public static String getName(String file, int level) {
		if (file == null)
			return null;
		// level can start at 1 for file name, 2 is for parent directory name
		String[] tokens = file.split(File.separator);
		if (tokens.length < level)
			throw new IllegalArgumentException(
					String.format(
							"Name for micrograph is taken from level %s, invalid path ",
							level, file));
		String name = tokens[tokens.length - level];
		if (level == 1) {
			int pos = name.lastIndexOf('.');
			if (pos != -1)
				name = name.substring(0, pos);

			pos = tokens[0].lastIndexOf('@');
			if (pos != -1) {
			    name = tokens[0].substring(0, pos) + "_at_" + name;
			}
		}
		return name;
	}

	public String getPosFile() {
		return posfile;
	}

	/* Load width and height after loaded ImagePlus */
	public void loadDimensions() {
                
		if (imp != null) {
			width = imp.getWidth();
			height = imp.getHeight();
		} else {
			try {
				ImageGeneric img = new ImageGeneric(file); // this reads the header
				width = img.getXDim();
				height = img.getYDim();
				img.destroy();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
        
        

	public ImagePlus getImagePlus() {
		try {

			if (imp == null) {
				imp = XmippImageConverter.loadImage(file);
				if (imp == null)
					imp = new ImagePlus(file);
				
			}
			return imp;
		} catch (Exception e) {
			ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public ImagePlus getImagePlus(List<IJCommand> filters) {
		try {
			if (filters.isEmpty())
				return getImagePlus();

			if (imp == null) {
				ImageGeneric ig = new ImageGeneric(file);
				ig.read(ImageGeneric.FIRST_IMAGE);

				// Smooth filter should be the first one
				// because it is applied in Xmipp
				for (IJCommand f : filters)
					if (f.getCommand().equals(ParticlePicker.xmippsmoothfilter)) {
						ig.convert2Datatype(ImageGeneric.UChar);
						ImageGeneric igsmooth = new ImageGeneric(
								ImageGeneric.UChar);
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
			// this filter was applied
			if (!f.getCommand().equals(ParticlePicker.xmippsmoothfilter)) 
			{
				IJ.run(imp, f.getCommand(), f.getOptions());
			}
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
        
    public int getWidth()
    {
        if(width == 0)
            loadDimensions();
        return width;
    }
    
    public int getHeigth()
    {
        if(height == 0)
            loadDimensions();
        return height;
    }

	public void releaseImage()
	{
		imp = null;
		
	}

    public List<? extends PickerParticle> getParticleList(){
        return null;
    };
}
