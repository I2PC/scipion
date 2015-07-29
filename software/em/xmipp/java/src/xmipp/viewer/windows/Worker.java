/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.windows;

import ij.ImagePlus;
import java.util.ArrayList;
import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MDRow;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;

/**
 *
 * @author airen
 */
public class Worker implements Runnable
	{
		public static final int STATS = 0;
		public static final int PCA = 1;
		public static final int FSC = 2;
		public String message;
		/** Constructor selecting operation */
		private int op; // store operation
		private MetaData imagesmd;
        private GalleryJFrame frame = null;

		public Worker(int operation, boolean selected, GalleryJFrame frame)
		{
            this.frame = frame;
                        
			op = operation;
            imagesmd = frame.data.getImagesMd(frame.gallery.getSelection(), selected);
			if (imagesmd.size() == 0)
				throw new IllegalArgumentException("No images available");
		}

		public void run()
		{
			try
			{
				switch (op)
				{
					case STATS:
						computeStatsImages();
						break;
					case PCA:
						pca();
						break;
					case FSC:
						fsc();
						break;
				}
			}
			catch (Exception e)
			{
				XmippWindowUtil.releaseGUI(frame.getRootPane());
				XmippDialog.showException(frame, e);
				return;

			}
			XmippWindowUtil.releaseGUI(frame.getRootPane());
		}

		public String getMessage()
		{
			switch (op)
			{
			case STATS:
				return "Computing average and std images...";
			case PCA:
				return "Computing PCA...";
			case FSC:
				return "Computing FSC...";
			}
			return "";
		}
                
        private void computeStatsImages() throws Exception
        {
            ImageGeneric imgAvg = new ImageGeneric();
            ImageGeneric imgStd = new ImageGeneric();

            imagesmd.getStatsImages(imgAvg, imgStd, frame.data.useGeo(), frame.data.isWrap(), MDLabel.MDL_IMAGE);
            ImagePlus impAvg = XmippImageConverter.convertToImagePlus(imgAvg);
            ImagePlus impStd = XmippImageConverter.convertToImagePlus(imgStd);
            imgAvg.destroy();
            imgStd.destroy();

            XmippImageWindow winAvg = new XmippImageWindow(new ImagePlusLoader(impAvg), "AVG:", frame.data.parameters);
            XmippWindowUtil.setLocation(0.2f, 0.5f, winAvg, frame);
            winAvg.setVisible(true);
            XmippImageWindow winStd = new XmippImageWindow(new ImagePlusLoader(impStd), "STD:", frame.data.parameters);

            XmippWindowUtil.setLocation(0.8f, 0.5f, winStd, frame);
            winStd.setVisible(true);
        }

        public void pca() throws Exception
        {
            ImageGeneric image = new ImageGeneric();
            imagesmd.getPCAbasis(image, MDLabel.MDL_IMAGE);
            ImagePlus imp = XmippImageConverter.convertToImagePlus(image);
            imp.setTitle("PCA: " + frame.data.getFileName());
            ImagesWindowFactory.openXmippImageWindow(imp, frame.data.parameters);

        }

        public void fsc() throws Exception
        {
            FSCJFrame fscframe = new FSCJFrame(frame.data, imagesmd, MDLabel.MDL_IMAGE);
            XmippWindowUtil.centerWindows(fscframe, frame);
            fscframe.setVisible(true);
        }
        
       
        
        
        public MDRow[] getImagesMd(MetaData md, int idlabel) {
            
            MDRow mdRow;
            ArrayList<MDRow> imagesmd = new ArrayList<MDRow>();
            int index = 0;
            String imagepath;
            for (long id : md.findObjects()) {
                if (frame.data.isEnabled(index)) {
                    imagepath = md.getValueString(idlabel, id, true);
                    System.out.println(imagepath);
                    if (imagepath != null && ImageGeneric.exists(imagepath)) {
                        mdRow = new MDRow();
                        if (frame.data.useGeo()) 
                            md.getRow(mdRow, id);//copy geo info in mdRow
                        mdRow.setValueString(idlabel, imagepath);
                        imagesmd.add(mdRow);
                    }
                }
                index++;
            }
            return imagesmd.toArray(new MDRow[]{});
        }

	}


