
import xmipp.*;

public class HelloWorld {

    public static void main(String args[]) {
        try {
            String path = "/media/PENDRIVE/Ad5GLflagIIIa/a09250/down3_a09250_Periodogramavg.psd";

            ImageDouble image = new ImageDouble();
            //image.read(path);
            image.readPreview(path, 128, 128);
            image.printShape();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static void main_(String args[]) {
        try {
            // Loads image plus for volume.
            String volumeFile = "/home/juanjo/temp/mask.vol";
            System.out.println("volume File: " + volumeFile);

            ImageDouble image = new ImageDouble();
            image.read(volumeFile);

            double rot = 0;
            double tilt = 0;

            Projection p = new Projection();
            System.out.println("1. Reset.");
            p.reset(image.getWidth(), image.getHeight());
            /*
            System.out.println("Image  : (W=" + image.getWidth() + " H=" + image.getHeight() + ")");
            System.out.println("Project: (w=" + p.getWidth() + " h=" + p.getHeight() + ")");

            System.out.println("2. Project.");
            Projection.projectVolume(image, p,
            image.getHeight(), image.getWidth(),
            rot, tilt, 0);

            System.out.println("3. Set XmippOrigin.");
            p.setXmippOrigin();

            System.out.println("4. Write.");
            p.write("/home/juanjo/Desktop/projection.xmp");*/
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
