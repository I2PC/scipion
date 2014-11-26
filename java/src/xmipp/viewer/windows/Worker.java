/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.windows;

import ij.ImagePlus;
import xmipp.ij.commons.ImagePlusLoader;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.ij.commons.XmippImageWindow;
import xmipp.jni.ImageGeneric;
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

		public Worker(int operation, MetaData imagesmd, GalleryJFrame frame)
		{
			op = operation;
			this.imagesmd = imagesmd;
			if (imagesmd.findObjects().length == 0)
				throw new IllegalArgumentException("No images available");
                        this.frame = frame;
		}

		public void run()
		{
			try
			{
				switch (op)
				{
				case STATS:
					computeStatsImages(imagesmd);
					break;
				case PCA:
					pca(imagesmd);
					break;
				case FSC:
					fsc(imagesmd);
					break;
				}
				imagesmd.destroy();
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
                
                private void computeStatsImages(MetaData imagesmd) throws Exception
                {
                        ImageGeneric imgAvg = new ImageGeneric();
                        ImageGeneric imgStd = new ImageGeneric();

                        imagesmd.getStatsImages(imgAvg, imgStd, frame.data.useGeo(), frame.data.getRenderLabel());
                        ImagePlus impAvg = XmippImageConverter.convertToImagePlus(imgAvg);
                        ImagePlus impStd = XmippImageConverter.convertToImagePlus(imgStd);
                        imgAvg.destroy();
                        imgStd.destroy();

                        XmippImageWindow winAvg = new XmippImageWindow(new ImagePlusLoader(impAvg), "AVG: " + frame.data.getFileName());
                        XmippWindowUtil.setLocation(0.2f, 0.5f, winAvg, frame);
                        winAvg.setVisible(true);
                        XmippImageWindow winStd = new XmippImageWindow(new ImagePlusLoader(impStd), "STD: " + frame.data.getFileName());

                        XmippWindowUtil.setLocation(0.8f, 0.5f, winStd, frame);
                        winStd.setVisible(true);
                        imagesmd.destroy();
                }

                public void pca(MetaData imagesmd) throws Exception
                {
                        ImageGeneric image = new ImageGeneric();
                        imagesmd.getPCAbasis(image, frame.data.getRenderLabel());
                        ImagePlus imp = XmippImageConverter.convertToImagePlus(image);
                        imp.setTitle("PCA: " + frame.data.getFileName());
                        ImagesWindowFactory.openXmippImageWindow(frame, imp, false);
                        imagesmd.destroy();

                }

                public void fsc(MetaData imagesmd) throws Exception
                {
                        FSCJFrame fscframe = new FSCJFrame(frame.data, imagesmd);
                        XmippWindowUtil.centerWindows(fscframe, frame);
                        fscframe.setVisible(true);
                }

	}


