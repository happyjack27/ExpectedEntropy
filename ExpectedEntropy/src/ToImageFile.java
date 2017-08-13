import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.nio.*;
import javax.imageio.*;

public class ToImageFile {
	public Draws draws = null;
	
 	public void toPNG(String file, Draws fd, int img_width) {
 		draws = fd;

		try {
			double size = img_width;
			double fsaa = 1;
			double dpi = 1;

    		BufferedImage image =
	        		  new BufferedImage(
	        				  (int) (size), 
	        				  (int) (size), 
	        		          BufferedImage.TYPE_INT_ARGB
	        		          );
	        Graphics2D graphics1 = image.createGraphics();
	        

	        graphics1.setRenderingHint(RenderingHints.KEY_INTERPOLATION,RenderingHints.VALUE_INTERPOLATION_BICUBIC);
	        graphics1.setRenderingHint(RenderingHints.KEY_RENDERING,RenderingHints.VALUE_RENDER_QUALITY);
	        graphics1.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
	        graphics1.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
	        graphics1.setStroke(new BasicStroke((int)fsaa));
	        graphics1.setFont(new Font("Arial",0,20));
	        draws.draw(graphics1,(int)size,(int)size);

			
    		BufferedImage bufferedImage2 =
	        		  new BufferedImage(
	        				  (int)(size*2*dpi), 
	        				  (int)(size*2*dpi), 
	        		          BufferedImage.TYPE_INT_ARGB
	        		          );
	        Graphics2D graphics2 = bufferedImage2.createGraphics();
	        graphics2.setRenderingHint(RenderingHints.KEY_INTERPOLATION,RenderingHints.VALUE_INTERPOLATION_BICUBIC);
	        graphics2.setRenderingHint(RenderingHints.KEY_RENDERING,RenderingHints.VALUE_RENDER_QUALITY);
	        graphics2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
	        graphics2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
	        graphics2.drawImage(image
		    		,0,0,(int)(size*2*dpi),(int)(size*2*dpi)
		    		,null);


    		BufferedImage bufferedImage =
	        		  new BufferedImage(
	        				  (int)(size*dpi), 
	        				  (int)(size*dpi), 
	        		          BufferedImage.TYPE_INT_ARGB
	        		          );
	        Graphics2D graphics = bufferedImage.createGraphics();
	        graphics.setRenderingHint(RenderingHints.KEY_INTERPOLATION,RenderingHints.VALUE_INTERPOLATION_BICUBIC);
	        graphics.setRenderingHint(RenderingHints.KEY_RENDERING,RenderingHints.VALUE_RENDER_QUALITY);
	        graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
	        graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

	        graphics.drawImage(bufferedImage2
		    		,0,0,(int)(size*dpi),(int)(size*dpi)
		    		,null);

			System.out.println("writting "+file);
			//ImageIO.write(bufferedImage, "png", new File(file));
			ImageIO.write(image, "png", new File(file));
		} catch (Exception e) {
			System.out.println("e "+e);
			e.printStackTrace();
		}
 		
 	}

}
