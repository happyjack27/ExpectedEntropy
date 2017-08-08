import java.io.*;
import java.util.Vector;

import distributions.CategoricalDistribution;

public class ProcessArgs {
	static int OUTPUT_FILE = 0;
	static int RESOLUTION = 1;
	static int COUNTS = 2;
	static int MARGINAL_CATS = 3;
	
	public static void main(String[] args) {
		args = "out.csv 10000 162,8,103,5,14,32,77,231,11,154,122,4,63,7,11,53,65,352,12,4 2".split(" ");
		boolean show_expectations = true;
		boolean show_curves = true;
		boolean use_prior_for_sampling = false;
		//args = "out.csv 1000 5,3,2,5,12,3 2".split(" ");
		//args = "1000 162,8,103,5,14,32,77,231,11,154,122,4,63, 7, 11, 53, 65, 352, 12, 4},		
		if( args == null || args.length < 4) {
			System.out.println("Usage: ");
			System.out.println("  java -jar catsnbees.jar output_file resolution counts marginal_categories [show expectations (true/false)] [show curves (true/false)] ");
			System.out.println("  example: java -jar catsnbees.jar out.csv 1000 51,3,20,5,12,3 2");
			System.out.println("  will use a marginal distribution of {5+3+2,  5+12+3}");
			System.out.println("       and a marginal distribution of {5+5, 3+12, 2+3}");
			return;
		}
		if( args.length > 4) {
			show_expectations = args[4].toLowerCase().equals("true");
		}
		if( args.length > 5) {
			show_curves = args[5].toLowerCase().equals("true");
		}
		int resolution = Integer.parseInt(args[RESOLUTION]);
		int cats = Integer.parseInt(args[MARGINAL_CATS]);
		String[] cs = args[COUNTS].split(",");
		//System.out.print("counts: ");
		Integer[] counts = new Integer[cs.length];
		for( int i = 0; i < counts.length; i++) {
			counts[i] = Integer.parseInt(cs[i]);
			//System.out.print(""+counts[i]+", ");
		}
		//System.out.println();
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
		cat.use_prior_for_sampling = use_prior_for_sampling;
		
		Vector<double[]>[] entropyCurves = cat.getAllEntropiesVectors(
				counts, margin1, margin2, resolution
		); 
		
		String[] names = new String[]{
				"H(X;Y)",
				"H(X)",
				"H(Y)",
				"H(X|Y)",
				"H(Y|X)",
				"I(X;Y)",
		};
		
		//now cat out the file.
		try {
			File f = new File(args[OUTPUT_FILE]);
			FileOutputStream fis = new FileOutputStream(f);
			
			StringBuffer sb = new StringBuffer();
			for( int i = 0; i < entropyCurves.length; i++) {
				Vector<double[]> vs = entropyCurves[i];
				if( show_expectations) {
					double expectation = cat.integrateSortedEntropyCurve(vs);
					sb.append("\"E("+names[i]+")\", ");
					appendLine(sb,new double[]{1,expectation});
				}
				if( show_curves) {
					for( int j = 0; j < vs.size(); j++) {
						sb.append("\""+names[i]+"\", ");
						appendLine(sb,vs.get(j));
					}
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
