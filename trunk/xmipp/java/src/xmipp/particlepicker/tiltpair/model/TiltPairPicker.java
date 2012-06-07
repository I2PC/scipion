package xmipp.particlepicker.tiltpair.model;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import xmipp.particlepicker.Family;
import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePicker;
import xmipp.particlepicker.training.model.FamilyState;
import xmipp.particlepicker.training.model.TrainingMicrograph;
import xmipp.utils.XmippMessage;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import java.util.Hashtable;

public class TiltPairPicker extends ParticlePicker
{

	protected List<UntiltedMicrograph> micrographs;
	private Family family;
	
	

	public TiltPairPicker(String selfile, String outputdir, FamilyState state)
	{
		super(selfile, outputdir, state);
		this.micrographs = new ArrayList<UntiltedMicrograph>();
		family = families.get(0);
		loadData();
	}
	

	private void loadData()
	{
		try
		{
			MetaData md = new MetaData(selfile);
			// md.readPlain(pairsfile, "image tilted_image");
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
				loadMicrographData(untiltedmicrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", selfile));

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void loadMicrographData(UntiltedMicrograph micrograph)
	{
		try
		{
			int x, y;
			UntiltedParticle up;
			TiltedParticle tp;
			String filename = getOutputPath(micrograph.getPosFile());
			if (!new File(filename).exists())
				return;

			MetaData md = new MetaData(filename);
			for (long id : md.findObjects())
			{
				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				up = new UntiltedParticle(x, y, micrograph, family);
				micrograph.addParticle(up);
			}
			filename = getOutputPath(micrograph.getTiltedMicrograph().getPosFile());
			md = new MetaData(filename);
			int i = 0;
			long[] ids = md.findObjects();
			for (long id : ids)
			{
				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				up = micrograph.getParticles().get(i);
				tp = new TiltedParticle(x, y, up);
				up.setTiltedParticle(tp);
				micrograph.getTiltedMicrograph().addParticle(tp);
				i++;
			}
			micrograph.initAligner();

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public int getSize()
	{
		return family.getSize();
	}

	public Color getColor()
	{
		return family.getColor();
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

	public void setColor(Color color)
	{
		family.setColor(color);
		setChanged(true);
	}

	public void resetMicrograph(UntiltedMicrograph m)
	{
		m.reset();
		setChanged(true);
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
		long id;

		try
		{
			MetaData md, md2, anglesmd;
			TiltedParticle tp;
			anglesmd = new MetaData(selfile);

			Hashtable<String, Long> micrographsDict = new Hashtable<String, Long>();
			for (long mid : anglesmd.findObjects())
				micrographsDict.put(anglesmd.getValueString(MDLabel.MDL_MICROGRAPH, mid), mid);

			for (UntiltedMicrograph m : micrographs)
			{
				if (!m.hasData())
					new File(m.getPosFile()).delete();
				else
				{

					md = new MetaData();
					md2 = new MetaData();
					id = micrographsDict.get(m.getFile());
					anglesmd.setValueDouble(MDLabel.MDL_ANGLE_Y, (double) m.getUntiltedAngle(), id);
					anglesmd.setValueDouble(MDLabel.MDL_ANGLE_Y2, (double) m.getTiltedAngle(), id);
					anglesmd.setValueDouble(MDLabel.MDL_ANGLETILT, (double) m.getTiltAngle(), id);

					for (UntiltedParticle p : m.getParticles())
					{
						tp = p.getTiltedParticle();
						if (tp != null)
						{
							id = md.addObject();
							md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
							md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);

							id = md2.addObject();
							md2.setValueInt(MDLabel.MDL_XINT, tp.getX(), id);
							md2.setValueInt(MDLabel.MDL_YINT, tp.getY(), id);
						}
					}
					String template = family.getName() + "@%s";
					md.write(String.format(template, getOutputPath(m.getPosFile())));
					md2.write(String.format(template, getOutputPath(m.getTiltedMicrograph().getPosFile())));
					anglesmd.write(selfile);
				}
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void setSize(int size)
	{
		family.setSize(size);
		setChanged(true);
	}

	public Family getFamily()
	{
		return family;
	}

	@Override
	public int getManualParticlesNumber(Family f)
	{
		int count = 0;
		for (UntiltedMicrograph um : micrographs)
			count += um.getParticles().size();
		return count;
	}

	public void importData(UntiltedMicrograph um, String file)
	{
		try
		{
			int x, y;
			UntiltedParticle up;
			MetaData md = new MetaData();
			md.readPlain(file, "Xcoor Ycoor");
			um.getParticles().clear();
			long[] ids = md.findObjects();
			for (long id : ids)
			{
				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				up = new UntiltedParticle(x, y, um, family);
				um.addParticle(up);
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	public void importData(TiltedMicrograph tm, String file)
	{
		try
		{
			int x, y;
			UntiltedParticle up;
			TiltedParticle tp;
			MetaData md = new MetaData();
			md.readPlain(file, "Xcoor Ycoor");
			int i = 0;
			long[] ids = md.findObjects();
			for (long id : ids)
			{
				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				up = tm.getUntiltedMicrograph().getParticles().get(i);
				tp = new TiltedParticle(x, y, up);
				up.setTiltedParticle(tp);
				tm.addParticle(tp);
				i++;
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}

	@Override
	public void exportParticles(Family family, String path)
	{
		throw new UnsupportedOperationException();

	}

	@Override
	public void importParticlesXmipp30Project(Family family, String absolutePath)
	{
		throw new UnsupportedOperationException(XmippMessage.getNotImplementedYetMsg());

	}

	@Override
	public void importParticlesFromXmipp24Project(Family family, String projectdir)
	{
		String ufile, tfile;

		for (UntiltedMicrograph um : micrographs)
		{
			ufile = getPosFileFromXmipp24Project(projectdir, um.getName());
			tfile = getPosFileFromXmipp24Project(projectdir, um.getTiltedMicrograph().getName());
			importParticlesFromXmipp24Files(um, ufile, tfile, false);
		}
	}

	public void importParticlesFromXmipp24Files(UntiltedMicrograph um, String ufile, String tfile, boolean checkfiles)
	{
		try
		{
			MetaData md = new MetaData();
			long[] ids;
			int x, y;
			UntiltedParticle up;
			TiltedParticle tp;
			TiltedMicrograph tm;

			um.getParticles().clear();
			if (!new File(ufile).exists())
			{
				if(checkfiles)
					throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("file", ufile));
				else 
					return;
			}
			um.setPosFileFromXmipp24(ufile);
			md.readPlain(ufile, "Xcoor Ycoor");
			ids = md.findObjects();
			for (long id : ids)
			{
				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				up = new UntiltedParticle(x, y, um, family);
				um.addParticle(up);
			}

			tm = um.getTiltedMicrograph();
			tm.getParticles().clear();
			if (!new File(tfile).exists())
			{
				if(checkfiles)
					throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("file", tfile));
				else 
					return;
			}
			tm.setPosFileFromXmipp24(tfile);
			md.readPlain(tfile, "Xcoor Ycoor");
			int i = 0;
			ids = md.findObjects();
			for (long id : ids)
			{
				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				up = um.getParticles().get(i);
				tp = new TiltedParticle(x, y, up);
				up.setTiltedParticle(tp);
				tm.addParticle(tp);
				i++;
			}

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}
}
