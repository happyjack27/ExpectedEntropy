import java.util.Date;
import java.util.Vector;

import javax.swing.JOptionPane;

import distributions.*;
import distributions.interfaces.PosteriorDistribution;
import util.Functions;
import org.apache.commons.math3.util.*;

public class Main {
	public static int verbosity = 0;
	
	public static void log(int level, String s) {
		if( verbosity < level) {
			return;
		}
		System.out.println(s);
	}
	public static int[] getSample() {
		int a = FastMath.random() > 0.5 ? 1 : 0;
		int b = FastMath.random() > 0.5 ? 1 : 0;
		int c = FastMath.random() > 0.5 ? 1 : 0;
		int d = b * c;
		int e = b + c > 0 ? 1 : 0;
		int f = 0;
		return new int[]{a,b,c,d,e,f};
	}

	public static void main(String[] args) {
		int resolution = 102400;
		
		
		CategoricalDistribution cat = new CategoricalDistribution();
		cat.prior_is_on_H = true;
		cat.add_to_all_cats = 0;
		long m0 = new Date().getTime();
		
		//		return new double[]{total_h_ph,0-total_e,total_ph};

		
		
		double[] dd0;
		double[] dd1;
		double[] dd2;
		
		// new Integer{162, 8, 103, 5, 14, 32, 77, 231, 11, 154, 122, 4, 63, 7, 11, 53, 65, 352, 12, 4}
		// new int[][]{new int[]{0,1,2,3,4,5,6,7,8,9}, new int[]{10,11,12,13,14,15,16,17,18,19}}
		// 
		
		dd0 = cat.getSummaryStats(new Integer[]{797,693}, resolution);
		System.out.println(dd0[0]);
		//797 693
		cat.use_prior_for_sampling = false;
		double ce3 = cat.getConditionalEntropy(
				//new Integer[]{80,176,175,105,139,94,156,129,112,86,116},
				//new int[][]{new int[]{0,1,2,3,4,5}, new int[]{6,7,8,9,10,11}},
				new Integer[]{162, 8, 103, 5, 14, 32, 77, 231, 11, 154, 122, 4, 63, 7, 11, 53, 65, 352, 12, 4},
				new int[][]{new int[]{0,1,2,3,4,5,6,7,8,9}, new int[]{10,11,12,13,14,15,16,17,18,19}},
				/*
				new int[][]{
					new int[]{0,10},
					new int[]{1,11},
					new int[]{2,12},
					new int[]{3,13},
					new int[]{4,14},
					new int[]{5,15},
					new int[]{6,16},
					new int[]{7,17},
					new int[]{8,18},
					new int[]{9,19},
					},
					*/
				resolution);
		JOptionPane.showMessageDialog(null, ""+ce3);
		System.out.println(ce3);
		
		System.exit(0);
		
		for( int i = 0; i < 20; i++) {
			dd0 = cat.getSummaryStats(new Integer[]{i*3,i,i*2,i*2,i,i*3}, resolution); //joint counts
			dd1 = cat.getSummaryStats(new Integer[]{i*4,i*4,i*4}, resolution); //counts of answers on question x
			dd2 = cat.getSummaryStats(new Integer[]{i*6,i*6}, resolution); //counts on dem or rep
			double ce2 = cat.getConditionalEntropy(new Integer[]{i*3,i,i*2,i*2,i,i*3}, 
					new int[][]{new int[]{0,1},new int[]{2,3},new int[]{4,5}}
					, resolution);
			double mi = dd1[0]+dd2[0]-dd0[0];
			double ce = dd0[0]-dd1[0];
			System.out.println(i+" : "+dd0[0]+", "+dd1[0]+","+dd2[0]+" | "+ce+", "+ce2);
		}
		System.exit(0);
		
		for( int i = 0; i < 100; i++) {
			dd0 = cat.getSummaryStats(new Integer[]{i,i,0,0}, resolution);
			dd1 = cat.getSummaryStats(new Integer[]{i,i}, resolution);
			dd1[0] *= 2;
			for( int j = 0; j < dd0.length; j++) {
				System.out.print(dd0[j]+", "+dd1[j]+", "+(dd0[j]-dd1[j]));
			}
			System.out.println();
		}
		
		/*
		resolution = 1024;
		cat.setNumberOfCategories(4);
		int m = 2;
		Vector<double[]> pdf = cat.getEntropyPDF(new Integer[]{m,m,0,0}, resolution*resolution, resolution);
		Vector<double[]> pdf2 = cat.getEntropyPDF(new Integer[]{m,m}, resolution*resolution, resolution);
		long m1 = new Date().getTime();
		
		for( int i = 0; i < pdf.size(); i++) {
			double[] dd = pdf.get(i);
			double[] dd2 = pdf2.get(i);
			System.out.print((dd[0]/FastMath.log(2))+", "+dd[1] +", "+ (2.0*dd2[0]/FastMath.log(2))+", "+dd2[1]);
			for( int j = 0; j < dd.length; j++) {
				//System.out.print(dd[j]+", ");
			}
			System.out.println();
		}
		*/
		
		//System.out.println("ms: "+(m1-m0));
		System.exit(0);
		
		/*
		double d;
		d = Functions.getExpectedEntropy(new Integer[]{0,0}, 1, 1)/FastMath.log(2.0);
		System.out.println("d: "+d);
		d = Functions.getExpectedEntropy(new Integer[]{1,0}, 1, 1)/FastMath.log(2.0);
		System.out.println("d: "+d);
		d = Functions.getExpectedEntropy(new Integer[]{1,1}, 1, 1)/FastMath.log(2.0);
		System.out.println("d: "+d);
		d = Functions.getExpectedEntropy(new Integer[]{10,1}, 1, 1)/FastMath.log(2.0);
		System.out.println("d: "+d);
		d = Functions.getExpectedEntropy(new Integer[]{100,1}, 1, 1)/FastMath.log(2.0);
		System.out.println("d: "+d);
		System.exit(0);
		*/
		
		//there is the input - the results of the sample trials, and the bayesian prior. 
		int n = 0;
		double p = 0;
		
		boolean prior_is_on_H = true;
		
		BernoulliDistribution posterior = new BernoulliDistribution();
		posterior.prior_is_on_H = prior_is_on_H;

		
		int total = n;
		int zeros = (int)(p*(double)n);
		double[] dd = showEntropyDensityForBernoulli(zeros,total-zeros, posterior,resolution);
		System.out.println("e: "+dd[0]+" h:"+dd[1]);
	}
		
	public static double[] showEntropyDensityForBernoulli( int zeros, int ones, BernoulliDistribution posterior, int num_sample_points ) {
		//initialize our array to collect results at different points.
		double[][] results = new double[num_sample_points][4];
		double sump = 0;
		
		//now fill up our array with the results
		for( int i = 0; i < num_sample_points; i++) {
			double theta = ((double)i) / (num_sample_points*2);
			
			//this produces samples tat are more uniform across entropy
			double t2 = 1.0/FastMath.log(1-theta);
			theta = 1.0/FastMath.log(theta);
			theta = theta/(theta+t2);
			
			results[i] = posterior.getEntropyAndProbabilityAtTheta(theta, new Integer[]{zeros, ones});
			sump += results[i][1];
		}
		
		//and print them out
		double last_h = 0;
		double last_ph = 0;
		double total_h_ph = 0;
		double total_ph = 0;
		double total_e = 0;
		for( int i = 0; i < num_sample_points; i++) {
			double h = results[i][0]/FastMath.log(2);
			double ph = results[i][1]/sump;
			double delta_h = h-last_h;
			double delta_ph = ph-last_ph;
			double delta_h_ph = 0;
			double delta_p = 0;
			delta_p = delta_h*(last_ph+ph)/2.0;
			total_ph += delta_p;
			last_ph = ph;
			last_h = h;
		}
		last_h = 0;
		last_ph = 0;
		for( int i = 0; i < num_sample_points; i++) {
			double h = results[i][0]/FastMath.log(2);
			double ph = results[i][1]/(sump*total_ph);
			double delta_h = h-last_h;
			double delta_ph = ph-last_ph;
			double delta_h_ph = 0;
			double delta_p = 0;
			delta_p = delta_h*(last_ph+ph)/2.0;
			if( delta_h != 0) {
				double ph_h = delta_ph/delta_h;
				double ph0 = last_ph - ph_h*last_h;
				double int_0 = integral(ph0,ph_h,last_h);// ph_h*(last_h*last_h)/2.0 + ph0*last_h;
				double int_1 = integral(ph0,ph_h,h);//ph_h*(h*h)/2.0 + ph0*h;
				delta_h_ph = int_1 - int_0; 
				
				
				total_h_ph += delta_h_ph;
				
				if( ph_h != 0) {
					double delta_e = integral_e_of_e(ph0,ph_h,h) - integral_e_of_e(ph0,ph_h,last_h);
					if( delta_e != delta_e) {
						log(1,"nan: "+ph0+","+ph_h+","+h+","+last_h);
						log(1," : "+(ph0+ph_h*h));
						log(1," : "+(ph0+ph_h*last_h));
					} else {
						total_e -= delta_e;
					}
				}
	
			}
			last_ph = ph;
			last_h = h;
			
			//log(1,h+", "+ph);
			log(1,h+", "+ph);//+", "+delta_h_ph+", "+delta_p);
			
		}
		total_e /= FastMath.log(2);
		//log(1,"E(H) = "+(total_h_ph/total_ph));
		log(1,"E(H) = "+(total_h_ph));
		log(1,"H(H) = "+(0-total_e));
		return new double[]{total_h_ph,0-total_e};
	}
	public static double integral(double a, double b, double x) {
		return x*x*(2.0*b*x+3.0*a)/6.0;
	}
	
	public static double integral_e_of_e(double a, double b, double x) {
		//integral of (a+bx)log(a+bx)
		double y = a+b*x;
		if( y == 0) {
			y = Double.MIN_VALUE;
		}
		return y*y*(2.0*FastMath.log(y)-1.0) / (4.0*b);
	}
}
