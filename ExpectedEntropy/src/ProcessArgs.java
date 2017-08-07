import java.io.*;
import java.util.Vector;

import distributions.CategoricalDistribution;

public class ProcessArgs {
	static int OUTPUT_FILE = 0;
	static int RESOLUTION = 1;
	static int COUNTS = 2;
	static int MARGINAL_CATS = 3;
	
	public static void main(String[] args) {
		if( args == null || args.length < 4) {
			System.out.println("Usage: ");
			System.out.println("  java -jar bees.jar output_file resolution counts marginal_categories");
			System.out.println("  example: out.csv 1000 5,3,2,5,12,3 2");
			System.out.println("  will use a marginal distribution of {5+3+2,  5+12+3}");
			System.out.println("       and a marginal distribution of {5+5, 3+12, 2+3}");
			return;
		}
		int resolution = Integer.parseInt(args[RESOLUTION]);
		int cats = Integer.parseInt(args[MARGINAL_CATS]);
		String[] cs = args[COUNTS].split(",");
		Integer[] counts = new Integer[cs.length];
		for( int i = 0; i < counts.length; i++) {
			counts[i] = Integer.parseInt(cs[i]);
		}
		//now make the marginals
		int other_cats = counts.length / cats;
		int[][] margin1 = new int[cats][other_cats];
		int[][] margin2 = new int[other_cats][cats];
		{
			int c = 0;
			for( int i = 0; i < margin1.length; i++) {
				for( int j = 0; j < margin2.length; j++) {
					margin1[i][j] = c;
					margin2[j][i] = c;
					c++;
				}			
			}
		}
		
		
		CategoricalDistribution cat = new CategoricalDistribution();
		cat.prior_is_on_H = true;
		
		Vector<double[]>[] entropyCurves = null;
		double[] summaries = new double[entropyCurves.length];
		for( int i = 0; i < entropyCurves.length; i++) {
			summaries[i] = cat.integrateSortedEntropyCurve(entropyCurves[i]);
		}
		
		String[] names = new String[]{
				"H(X;Y)",
				"H(X)",
				"H(Y)",
				"H(X|Y)",
				"H(Y|X)",
				"I(X;Y)",
				
		};
		
		//Vector<double[]> vs = new Vector<double[]>();
		
		//now cat out the file.
		try {
			File f = new File(args[OUTPUT_FILE]);
			FileOutputStream fis = new FileOutputStream(f);
			
			StringBuffer sb = new StringBuffer();
			sb.append("\"EXPECTATIONS\", ");
			appendLine(sb,summaries);
			for( int i = 0; i < entropyCurves.length; i++) {
				Vector<double[]> vs = entropyCurves[i]; 
				for( int j = 0; j < vs.size(); j++) {
					sb.append("\""+names[i]+"\", ");
					appendLine(sb,vs.get(i));
				}
			}
			fis.write(sb.toString().getBytes());
			fis.flush();
			fis.close();
			System.out.println("Done.");
			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	public static void appendLine(StringBuffer sb, double[] dd) {
		for( int j = 0; j < dd.length; j++) {
			if( j > 0) {
				sb.append(", ");
			}
			sb.append(""+dd[j]);
		}
		sb.append("\n");
	}
}
