package xmipp.viewer.particlepicker.tiltpair.model;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;

import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.XmippMessage;
import xmipp.viewer.particlepicker.Format;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.training.model.Mode;


/**
 * Business object for Tilt Pair Picker GUI. Inherits from ParticlePicker 
 * @author airen
 *
 */
public class TiltPairPicker extends ParticlePicker
{

	protected List<UntiltedMicrograph> micrographs;
	private UntiltedMicrograph micrograph;

	public TiltPairPicker(String selfile, String outputdir, Mode state)
	{
		super(selfile, outputdir, state);

		for (UntiltedMicrograph um : micrographs)
			loadMicrographParticles(um);
	}

	public void loadData()
	{
		try
		{
			loadEmptyMicrographs();
			for (UntiltedMicrograph um : micrographs)
				loadMicrographParticles(um);

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void loadEmptyMicrographs()
	{
		try
		{
			MetaData md = new MetaData(selfile);
			// md.readPlain(pairsfile, "image tilted_image");
			if (micrographs == null)
				this.micrographs = new ArrayList<UntiltedMicrograph>();
			else
				micrographs.clear();
			UntiltedMicrograph untiltedmicrograph;
			TiltedMicrograph tiltedmicrograph;
			String image, tiltedimage;

			long[] ids = md.findObjects();
			for (long id : ids)
			{
				image = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
				tiltedimage = md.getValueString(MDLabel.MDL_MICROGRAPH_TILTED, id);
				tiltedmicrograph = new TiltedMicrograph(tiltedimage);
				untiltedmicrograph = new UntiltedMicrograph(image, tiltedmicrograph);
				tiltedmicrograph.setUntiltedMicrograph(untiltedmicrograph);
				micrographs.add(untiltedmicrograph);

			}
			md.destroy();
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", selfile));

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void loadMicrographParticles(UntiltedMicrograph micrograph)
	{
		String ufile = getOutputPath(micrograph.getPosFile());
		String tfile = getOutputPath(micrograph.getTiltedMicrograph().getPosFile());
		if (Filename.exists(ufile) && Filename.exists(tfile))
			loadMicrographData(micrograph, ufile, tfile);
	}

	/* Return number of particles loaded */
	public String loadMicrographParticles(UntiltedMicrograph um, MetaData uMd, MetaData tMd)
	{
		String result = "";
		um.reset(this);
		UntiltedParticle up;
		TiltedParticle tp;
		TiltedMicrograph tm = um.getTiltedMicrograph();
		int x, y;
		long[] uIds = uMd.findObjects();
		long[] tIds = tMd.findObjects();
		long id;
		for (int i = 0; i < uIds.length; ++i)
		{
			// Add untilted particle
			id = uIds[i];
			x = uMd.getValueInt(MDLabel.MDL_XCOOR, id);
			y = uMd.getValueInt(MDLabel.MDL_YCOOR, id);
			up = new UntiltedParticle(x, y, um, this);
			um.addParticle(up);
			// Set tilted pair particle
			if (i < tIds.length)
			{
				id = tIds[i];
				x = tMd.getValueInt(MDLabel.MDL_XCOOR, id);
				y = tMd.getValueInt(MDLabel.MDL_YCOOR, id);
				if (x <= 0 || y <= 0)
				{
					result += String.format("Tilted particle at %s centered on %s,%s with negative coordinates dismissed.\n", getMicrograph()
							.getName(), x, y);
					continue;
				}
				tp = new TiltedParticle(x, y, up);
				up.setTiltedParticle(tp);
				tm.addParticle(tp);
			}
			else
				result += String.format("Particle at %s centered on %s,%s without tilted pair.\n", getMicrograph().getName(), x, y);
		}

		um.initAligner();

		return result;
	}// loadMicrographs

	public void loadMicrographData(UntiltedMicrograph um, String uPosFile, String tPosFile)
	{
		try
		{
			MetaData uMd = new MetaData(getParticlesBlock(uPosFile));
			MetaData tMd = new MetaData(getParticlesBlock(tPosFile));
			loadMicrographParticles(um, uMd, tMd);
			uMd.destroy();
			tMd.destroy();
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public int getNextFreeMicrograph()
	{
		int count = 0;
		for (UntiltedMicrograph m : micrographs)
		{
			if (m.hasData())
				return count;
			count++;
		}
		return -1;
	}

	public List<UntiltedMicrograph> getMicrographs()
	{
		return micrographs;
	}

	public void resetMicrograph(UntiltedMicrograph m)
	{
		m.reset(this);
		saveData(m);
	}

	public int getUntiltedNumber()
	{
		int count = 0;
		for (UntiltedMicrograph um : micrographs)
			count += um.getParticles().size();
		return count;
	}

	public int getTiltedNumber()
	{
		int count = 0;
		for (UntiltedMicrograph um : micrographs)
			count += um.getTiltedMicrograph().getParticles().size();
		return count;
	}

	public void saveData()
	{
		super.saveData();
		saveMicrographAngles(micrograph);
	}
	
	private void saveMicrographAngles(UntiltedMicrograph m)
	{
		try
		{
			MetaData anglesmd;
			anglesmd = new MetaData(selfile);
			
			long micId = -1; //Current micrograph id
			long [] ids = anglesmd.findObjects();
			String micFile = micrograph.getFile();
			
			for (long mid: ids)
			{
				if (micFile.equals(anglesmd.getValueString(MDLabel.MDL_MICROGRAPH, mid)))
				{
					micId = mid;
					break;
				}
			}	
			
			if (micId == -1)
				throw new Exception("Micrograph " + micrograph.getFile() + " was not found in metadata");
			
			anglesmd.setValueDouble(MDLabel.MDL_ANGLE_Y,  m.getUntiltedAngle(), micId);
			anglesmd.setValueDouble(MDLabel.MDL_ANGLE_Y2, m.getTiltedAngle(), micId);
			anglesmd.setValueDouble(MDLabel.MDL_ANGLE_TILT, m.getTiltAngle(), micId);
			
			anglesmd.writeBlock(selfile);
			saveData(m);
		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
		
	}

	@Override
	public void saveData(Micrograph m)
	{
		UntiltedMicrograph um = (UntiltedMicrograph) m;
		String file = getOutputPath(um.getPosFile());
		long id;

		try
		{

			if (!m.hasData()){
				new File(file).delete();
			}
			else
			{
				TiltedParticle tp;

				MetaData mdU = new MetaData(); // untilted micrograph particles
				MetaData mdT = new MetaData(); // tilted micrograph particles

				for (UntiltedParticle p : um.getParticles())
				{
					tp = p.getTiltedParticle();
					if (tp != null)
					{
						id = mdU.addObject();
						mdU.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
						mdU.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);

						id = mdT.addObject();
						mdT.setValueInt(MDLabel.MDL_XCOOR, tp.getX(), id);
						mdT.setValueInt(MDLabel.MDL_YCOOR, tp.getY(), id);
					}
				}

				mdU.write(getParticlesBlock(file));
				file = getOutputPath(um.getTiltedMicrograph().getPosFile());
				mdT.write(getParticlesBlock(file));
				
				mdU.destroy();
				mdT.destroy();
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public String importParticlesFromFolder(String path, Format f, float scale, boolean invertx, boolean inverty)
	{
		if (f == Format.Auto)
			f = detectFormat(path);
		if (f == Format.Unknown)
			throw new IllegalArgumentException("Unable to detect format");

		String uFn, tFn;
		String result = "";

		for (UntiltedMicrograph um : micrographs)
		{
			uFn = getImportMicrographName(path, um.getFile(), f);
			tFn = getImportMicrographName(path, um.getTiltedMicrograph().getFile(), f);
			if (Filename.exists(uFn) && Filename.exists(tFn))
			{
				result += importParticlesFromFiles(uFn, tFn, f, um, scale, invertx, inverty);
				saveData(um);
			}
		}

		return result;
	}// function importParticlesFromFolder

	public String importParticlesFromFiles(String uPath, String tPath, Format f, UntiltedMicrograph um, float scale, boolean invertx, boolean inverty)
	{
		MetaData uMd = new MetaData();
		fillParticlesMdFromFile(uPath, f, um, uMd, scale, invertx, inverty);
		MetaData tMd = new MetaData();
		fillParticlesMdFromFile(tPath, f, um.getTiltedMicrograph(), tMd, scale, invertx, inverty);

		String result = loadMicrographParticles(um, uMd, tMd);
		uMd.destroy();
		tMd.destroy();
		return result;
	}// function importParticlesFromFiles

	public String getImportMicrographName(String path, String filename, Format f)
	{
		String base = Filename.removeExtension(Filename.getBaseName(filename));
		switch (f)
		{
		case Xmipp24:
			return Filename.join(path, base, base + ".raw.Common.pos");
		case Xmipp24a:
			return Filename.join(path, base + ".raw.Common.pos");
		case Xmipp24b:
			return Filename.join(path, base.replaceAll("down[0-9]_", ""), base + ".raw.Common.pos");
		case Xmipp24c:
			return Filename.join(path, base + ".raw.pos");
		case Xmipp30:
			return Filename.join(path, base + ".pos");
		case Xmipp301:
			return Filename.join(path, base + ".pos");
		case Eman:
			return Filename.join(path, base + "_ptcls.box");

		default:
			return null;
		}
	}

	public UntiltedMicrograph getMicrograph()
	{
		return micrograph;
	}

	@Override
	public void setMicrograph(Micrograph m)
	{
		this.micrograph = (UntiltedMicrograph) m;

	}

	@Override
	public boolean isValidSize(int size)
	{
		UntiltedMicrograph um = getMicrograph();
		for (UntiltedParticle p : um.getParticles())
			if (!getMicrograph().fits(p.getX(), p.getY(), size))
				return false;
		for (TiltedParticle p : um.getTiltedMicrograph().getParticles())
			if (!um.getTiltedMicrograph().fits(p.getX(), p.getY(), size))
				return false;
		return true;
	}

}// class TiltPairPicker
