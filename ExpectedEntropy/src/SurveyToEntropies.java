import java.io.*;
import java.util.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.IntegerSequence;

import distributions.CategoricalDistribution;

public class SurveyToEntropies {
	public static final int OUTPUT_FILE = 0;
	public static final int INPUT_FILE = 1;
	public static final int VARS = 2;

	public static boolean use_bayesian = true;
	public static boolean full_cats = false;
	public static boolean prior_is_on_H = false;

	public static boolean cap_at_zero = true;
	
	public static boolean permuted = false;
	
	public static int bayesian_resolution = 1024;
	
	public static String delimiter = "\t";
	public static int numVars = 10;
	
	public static int[] choice_numbers = null;
	
	public static void main(String[] args) {
		try {
			if( args.length == 0) {
				args = new String[]{
						"",
						"data/first10.txt",
						""+5
				};
			}
			numVars = Integer.parseInt(args[VARS]);
			File f = new File(args[INPUT_FILE]);
			
			FileInputStream fis = new FileInputStream(f);
			BufferedReader buf = new BufferedReader(new InputStreamReader(fis));
			Vector<int[]> answers = new Vector<int[]>();
			HashMap<String,Integer>[] choices = new HashMap[numVars];
			for( int i = 0; i < numVars; i++) {
				choices[i] = new HashMap<String,Integer>();
			}
			choice_numbers = new int[numVars];
			
			
			System.out.println("reading file...");
			Thread.sleep(1);
			while(buf.ready()) {
				String[] ss = buf.readLine().split(delimiter);
				int[] answer = new int[numVars];
				for( int i = 0; i < numVars; i++) {
					Integer n = choices[i].get(ss[i]);
					if( n == null) {
						n = choice_numbers[i]++;
						choices[i].put(ss[i],n);
					}
					answer[i] = n;
				}
				answers.add(answer);
				Thread.sleep(1);
			}
			System.out.println("read "+answers.size()+" answers");
			
			System.out.println("computing entropies...");
			double[] unions = null;
			if( use_bayesian) {
				unions = computeAllJointEntropiesBayeisan(answers);
			} else {
				unions = computeAllJointEntropies(answers);
			}
			System.out.println("converting...");
			
			double[] intersections = unionsToIntersections(unions);
			for( int i = 0; i < intersections.length; i++) {
				if( Integer.bitCount(i) > 2) {
					//intersections[i] = 0;
				}
			}			
			unions = unionsToIntersections(intersections);
			intersections = unionsToIntersections(unions);
			
			intersections[0] = unions[unions.length-1];
			
			System.out.println("intersections");
			for( int i = 0; i < intersections.length; i++) {
				System.out.println(intersections[i]+",");
			}
			
			System.out.println("unions");
			for( int i = 0; i < unions.length; i++) {
				//System.out.println(unions[i]);
			}
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public static double[] computeAllJointEntropies(Vector<int[]> answers) {
		int size = 0x01 << numVars;
		double[] ret = new double[size];
		ret[0] = 0;
		for( int i = 1; i < size; i++) {
			ret[i] = computeEntropyMLE(createJointBin(answers,i));
		}
		return ret;
	}
	public static double[] computeAllJointEntropiesBayeisan(Vector<int[]> answers) {
		CategoricalDistribution cat  = new CategoricalDistribution();
		cat.prior_is_on_H = prior_is_on_H;
		int size = 0x01 << numVars;
		double[] ret = new double[size];
		ret[0] = 0;
		for( int i = 1; i < size; i++) {
			Integer[] ii = createJointBin(answers,i);
			Vector<double[]> curve = cat.getEntropyCurveLogarithmic(ii,bayesian_resolution,permuted);
			//Vector<double[]> curve = cat.getEntropyCurveMultiplicative(ii,bayesian_resolution);
			ret[i] = cat.integrateSortedEntropyCurve(curve);
			if( Double.isInfinite(ret[i]) || Double.isNaN(ret[i]) || curve.size() == 0) {
				ret[i] = computeEntropyMLE(ii);
			}
			if( Double.isInfinite(ret[i]) || Double.isNaN(ret[i])) {
				ret[i] = 0;
			}
		}
		return ret;
	}
	
	//http://www.elegantcoding.com/2012/10/the-ubiquitous-patterns-of-pascals.html
	public static double[] unionsToIntersections(double[] in) {
		int size = 0x01 << numVars;
		double[] ret = new double[size];
		for( int i = 0; i < size; i++) {
			ret[i] = 0;
		}
		for( int i = 1; i < size; i++) {
			for( int j = 1; j < size; j++) {
				if( (i & j) == j) { //j is a subset of i
					ret[i] += (Integer.bitCount(j) % 2 == 0) ? -in[j] : in[j];
				}
			}
		}
		if( cap_at_zero) {
			for( int i = 1; i < size; i++) {
				if( ret[i] < 0) {
					ret[i] = 0;
				}
			}
		}
		return ret;
	}
	
	public static void populateHashMap(HashMap<String,Integer> bins, int bits, int init_shift, String s) {
		for( int shift = init_shift; shift < 16 && shift < numVars; shift++) {
			if( (bits & (0x01 << shift)) != 0) {
				for( int i = 0; i < choice_numbers[shift]; i++) {
					String s2 = s+i+",";
					populateHashMap( bins, bits, shift+1, s2); 
				}
				//s+= answer[shift] + ",";
				return;
			}
		}
		bins.put(s, 0);
	}
	
	public static Integer[] createJointBin(Vector<int[]> answers, int bits) {
		HashMap<String,Integer> bins = new HashMap<String,Integer>();
		
		if( use_bayesian && full_cats) {
			populateHashMap(bins,bits,0,"");
		}
		

		for( int[] answer : answers) {
			String s = "";
			for( int shift = 0; shift < 16 && shift < numVars; shift++) {
				if( (bits & (0x01 << shift)) != 0) {
					s+= answer[shift] + ",";
				}
			}
			Integer n = bins.get(s);
			if( n == null) {
				n = 0;
			}
			bins.put(s, n+1);
		}
		
		Integer[] ret = new Integer[bins.size()];
		int r = 0;
		for( Integer n : bins.values()) {
			ret[r++] = n;
		}
		return ret;
	}
	public static double computeEntropyMLE(Integer[] bins) {
		double[] ps = new double[bins.length];
		double tot = 0;
		for( int i = 0; i < bins.length; i++) {
			ps[i] = bins[i];
			tot += bins[i];
		}
		tot = 1.0/tot;
		double ret = 0;
		for( int i = 0; i < bins.length; i++) {
			ps[i] *= tot;
			ret -= ps[i]*FastMath.log(ps[i]);
		}
		return ret;
	}
}
