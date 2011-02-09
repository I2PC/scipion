
import xmipp.*;

public class ImageTest {

    public static void main(String args[]) {
try{
        String file = "/home/juanjo/temp/avg0.spi";

        ImageDouble image = new ImageDouble();
        image.readPreview(file, 80, 80, 0, 0);
        image.printShape();
}catch(Exception ex){
ex.printStackTrace();
}
    }
}
