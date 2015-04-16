package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import java.awt.Image;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.ServerSocket;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import xmipp.jni.ImageGeneric;


public class XmippUtil {

	private static XmippImageJ xij;

	public static XmippImageJ showImageJ(Tool tool) {
		if (IJ.getInstance() == null) {
		try {
                        xij = new XmippImageJ();

//			
                File tempFile = File.createTempFile("macros", ".txt");
                BufferedWriter writer = new BufferedWriter(new FileWriter(tempFile));
                writer.write("macro \"Particle Picker Tool - C0a0L18f8L818f\" {   }\nmacro \"Xmipp Micrograph Viewer Tool - C0a0L18f8L818f\" {  }");
                writer.close();
                IJ.run("Install...", "install=" + tempFile.getAbsolutePath());
                //Tool cannot be set if it is not a default one
             } catch (Exception ex) {
                Logger.getLogger(XmippUtil.class.getName()).log(Level.SEVERE, null, ex);
                throw new IllegalArgumentException(ex);
            }
		} else if (!xij.isVisible())
			xij.setVisible(true);
		return xij;
	}

	public static XmippImageJ getXmippImageJ() {
		return xij;
	}
        
        

	public static ImagePlus getImagePlus(String file) {
		try {

			ImageGeneric ig = new ImageGeneric(file);
			ig.read(ImageGeneric.FIRST_IMAGE);
			ImagePlus imp = XmippImageConverter.convertToImagePlus(ig);
			ig.destroy();

			return imp;
		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	public static Icon getImageIcon(Image imp, int width, int height)
	{

		Image image = imp.getScaledInstance(width, height, Image.SCALE_SMOOTH);
		Icon icon = new ImageIcon(image);

		return icon;
	}
        


    


    public static void copyFile(String source, String dest) throws IOException
    {
        copyFile(new File(source), new File(dest));
    }
    
    public static void copyFile(File source, File dest)
            throws IOException {
        InputStream input = null;
        OutputStream output = null;
        try {
            input = new FileInputStream(source);
            output = new FileOutputStream(dest);
            byte[] buf = new byte[1024];
            int bytesRead;
            while ((bytesRead = input.read(buf)) > 0) {
                output.write(buf, 0, bytesRead);
            }
        } finally {
            if(input != null)
                input.close();
            if(output != null)
                output.close();
        }
    }
    
    public static boolean isInPath(String executable)
    {
        String[] paths = System.getenv("PATH").split(Pattern.quote(File.pathSeparator));
        String exePath;
        File file, exeFile;
        for(String path: paths)
        {
            file = new File(path);
            if(file.isDirectory())
            {
                exePath = path + File.separator + executable;
                exeFile = new File(exePath);
                if(exeFile.exists())
                    return true;
            }
            else if(file.getName().equals(executable))
                return true;
                
        }
        return false;
    }
    
     public static String formatNumbers(String str)
    {
        Pattern p = Pattern.compile("(-?(\\d)+(\\.)?(\\d)*)");
        Matcher m = p.matcher(str);
        StringBuffer sb = new StringBuffer(str.length());
        while(m.find())
        {
            String number = m.group(1);
            m.appendReplacement(sb, String.format("%.2f", Double.parseDouble(number)));
        }
        m.appendTail(sb);
        return sb.toString();
    }
     
    public static int findFreePort() 
             throws IOException {
            ServerSocket server = new ServerSocket(0);
            int port = server.getLocalPort();
            server.close();
            return port;
    }
       
}
