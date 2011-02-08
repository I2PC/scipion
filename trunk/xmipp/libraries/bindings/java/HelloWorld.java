
import xmipp.*;

class HelloWorld {
public static void main(String args[]) {
        String path = "/home/juanjo/Desktop/imgs_Roberto/kk.mrcs";
        ImageDouble image = new ImageDouble();

        try {
            image.read(path);
            System.err.println(" *** Nimages: " + image.getNimages());
        } catch (Exception ex) {
            ex.printStackTrace();
            System.err.println(" *** Exception: " + ex.getMessage());
        }
    }
}
