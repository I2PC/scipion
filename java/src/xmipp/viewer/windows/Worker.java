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
import xmipp.jni.MDLabel;
import xmipp.jni.MDRow;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippMessage;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.scipion.ScipionGalleryData;
import xmipp.viewer.scipion.ScipionGalleryJFrame;

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
                private int renderLabel;

		public Worker(int operation, MetaData md, GalleryJFrame frame)
		{
                        this.frame = frame;
                        
			op = operation;
                        renderLabel = frame.data.getRenderLabel();
                        if(frame.data instanceof ScipionGalleryData)
                            renderLabel = getFirstXmippRenderLabel(md);
                        

                        imagesmd = getImagesMd(md, renderLabel);
                        imagesmd.print();
			if (imagesmd.findObjects().length == 0)
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
                
                private void computeStatsImages() throws Exception
                {
                        ImageGeneric imgAvg = new ImageGeneric();
                        ImageGeneric imgStd = new ImageGeneric();

                        imagesmd.getStatsImages(imgAvg, imgStd, frame.data.useGeo(), renderLabel);
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

                public void pca() throws Exception
                {
                        ImageGeneric image = new ImageGeneric();
                        imagesmd.getPCAbasis(image, renderLabel);
                        ImagePlus imp = XmippImageConverter.convertToImagePlus(image);
                        imp.setTitle("PCA: " + frame.data.getFileName());
                        ImagesWindowFactory.openXmippImageWindow(frame, imp, false);
                        imagesmd.destroy();

                }

                public void fsc() throws Exception
                {
                        FSCJFrame fscframe = new FSCJFrame(frame.data, imagesmd);
                        XmippWindowUtil.centerWindows(fscframe, frame);
                        fscframe.setVisible(true);
                }
                
                public int getFirstXmippRenderLabel(MetaData xmippmd) {
                    
                    int[] labels = xmippmd.getActiveLabels();
                    for(int i = 0; i < labels.length; i ++)
                        if(MetaData.isImage(labels[i]))
                        {
                            return labels[i];
                        }
                    throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("Xmipp render label"));
                }
                
                
                public MetaData getImagesMd(MetaData md, int idlabel) {
                    
                    MDRow mdRow = new MDRow();
                    MetaData imagesmd = new MetaData();
                    int index = 0;
                    String imagepath;
                    long id2;
                    for (long id : md.findObjects()) {
                        if (frame.data.isEnabled(index)) {
                            imagepath = md.getValueString(idlabel, id, true);
                            System.out.println(imagepath);
                            if (imagepath != null && ImageGeneric.exists(imagepath)) {
                                id2 = imagesmd.addObject();
                                if (frame.data.useGeo()) {
                                    md.getRow(mdRow, id);
                                    mdRow.setValueString(idlabel, imagepath);
                                    imagesmd.setRow(mdRow, id2);
                                } else {
                                    imagesmd.setValueString(idlabel, imagepath, id2);
                                }
                            }
                        }
                        index++;
                    }
                    mdRow.destroy();
                    return imagesmd;
                }

	}


