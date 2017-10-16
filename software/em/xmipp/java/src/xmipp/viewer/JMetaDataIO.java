/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.logging.Level;

import xmipp.ij.commons.IJCommand;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import static xmipp.viewer.particlepicker.ParticlePicker.getLogger;
import xmipp.viewer.particlepicker.training.model.AutomaticParticle;
import xmipp.viewer.particlepicker.training.model.ManualParticle;
import xmipp.viewer.particlepicker.training.model.SupervisedParticlePicker;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;

/**
 *
 * @author airen
 */
public class JMetaDataIO {
    
    
    public static void saveData(SupervisedPickerMicrograph m, String path)
	{
                //System.out.println("saving data on JMetaDataIO");
		try
		{
			if (!m.hasData())
				new File(path).delete();
			else
			{
                             FileWriter fstream = new FileWriter(path);
                             BufferedWriter out = new BufferedWriter(fstream);
                             out.write("# XMIPP_STAR_1 * \n#\ndata_header\nloop_\n");
                             String cformat = "_%s\n";
                             out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_PICKING_MICROGRAPH_STATE)));
                             out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_PICKING_AUTOPICKPERCENT)));
                             out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_COST)));
                             String vformat = "%20s";
                             out.write(String.format(vformat, m.getState().toString()));
                             out.write(String.format(vformat, m.getAutopickpercent()));
                             out.write(String.format(vformat, m.getThreshold()));
                             out.newLine();
                             out.write("data_particles\nloop_\n");
                             out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_XCOOR)));
                             out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_YCOOR)));
                             
                            for (ManualParticle p : m.getManualParticles())
                            {
                                    out.write(String.format(vformat, p.getX()));
                                    out.write(String.format(vformat, p.getY()));
                                    out.newLine();
                            }
                            if(m.hasAutomaticParticles())
                            {
                                out.write("data_particles_auto\nloop_\n");
                                out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_XCOOR)));
                                out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_YCOOR)));
                                out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_COST)));
                                out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_ENABLED)));
                                int enabled;
                                for (AutomaticParticle p : m.getAutomaticParticles())
                                {
                                        out.write(String.format(vformat, p.getX()));
                                        out.write(String.format(vformat, p.getY()));
                                        out.write(String.format(vformat, p.getCost()));
//                                        enabled = !p.isDeleted()? 1 : -1;
                                        enabled = !p.isUnavailable() ? 1 : -1;
                                        out.write(String.format(vformat, enabled));
                                        out.newLine();
                                }
                            }
                            out.close();
			}                        
            }
            catch (Exception e)
            {
                    getLogger().log(Level.SEVERE, e.getMessage(), e);
                    throw new IllegalArgumentException(e.getMessage());
            }

	}
    
       
    public static void saveConfig(SupervisedParticlePicker picker, String path)
	{
                //System.out.println("saving data on JMetaDataIO");
		long id;
		try
		{
			
                    FileWriter fstream = new FileWriter(path);
                    BufferedWriter out = new BufferedWriter(fstream);
                    out.write("# XMIPP_STAR_1 * \n#\ndata_properties\nloop_\n");
                    String cformat = "_%s\n";
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_MICROGRAPH)));
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_COLOR)));
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_PICKING_PARTICLE_SIZE)));
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_PICKING_AUTOPICKPERCENT)));
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_PICKING_TEMPLATES)));
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_PICKING_STATE)));
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_PICKING_MANUALPARTICLES_SIZE)));
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_PICKING_AUTOPARTICLES_SIZE)));
                    
                    String vformat = "%20s";
                    out.write(String.format(vformat, picker.getMicrograph().getName()));
                    out.write(String.format(vformat, picker.getColor().getRGB()));
                    out.write(String.format(vformat, picker.getSize()));
                    out.write(String.format(vformat, picker.getAutopickpercent()));
                    out.write(String.format(vformat, picker.getTemplatesNumber()));
                    out.write(String.format(vformat, picker.getMode().toString()));
                    out.write(String.format(vformat, picker.getManualParticlesNumber()));
                    out.write(String.format(vformat, picker.getAutomaticParticlesNumber()));
                    out.newLine();
                    String optionsFormat = "%80s";
                    out.write("data_filters\nloop_\n");
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_MACRO_CMD)));
                    out.write(String.format(cformat, MetaData.getLabelName(MDLabel.MDL_MACRO_CMD_ARGS)));
                    String options;
                    for (IJCommand f : picker.getFilters()) {
                        out.write(String.format(vformat, f.getMdCommand()));
                        options = f.getMdOptions();
                        out.write(String.format(optionsFormat, options));
                        out.newLine();
                        
                    }
                    out.close();
            }
            catch (Exception e)
            {
                    getLogger().log(Level.SEVERE, e.getMessage(), e);
                    throw new IllegalArgumentException(e.getMessage());
            }

	}
    
}
