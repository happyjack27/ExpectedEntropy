import java.util.Vector;

import distributions.*;
import distributions.interfaces.PosteriorDistribution;

public class Main {

	public static void main(String[] args) {
		
		//there is the input - the results of the sample trials, and the bayesian prior. 
		int total = 16;
		int zeros = 0;
		boolean prior_is_on_H = true;
		int num_sample_points = 256;
		
		BernoulliDistribution posterior = new BernoulliDistribution();
		posterior.prior_is_on_H = prior_is_on_H;

		
		showEntropyDensityForBernoulli(zeros,total-zeros, posterior,num_sample_points);
	}
		
	public static void showEntropyDensityForBernoulli( int zeros, int ones, BernoulliDistribution posterior, int num_sample_points ) {
		//initialize our array to collect results at different points.
		double[][] results = new double[num_sample_points][4];
		double sump = 0;
		
		//now fill up our array with the results
		for( int i = 0; i < num_sample_points; i++) {
			double theta = ((double)i) / (num_sample_points*2);
			
			//this produces samples tat are more uniform across entropy
			double t2 = 1.0/Math.log(1-theta);
			theta = 1.0/Math.log(theta);
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
			double h = results[i][0]/Math.log(2);
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
			double h = results[i][0]/Math.log(2);
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
						System.out.println("nan: "+ph0+","+ph_h+","+h+","+last_h);
						System.out.println(" : "+(ph0+ph_h*h));
						System.out.println(" : "+(ph0+ph_h*last_h));
					} else {
						total_e -= delta_e;
					}
				}
	
			}
			last_ph = ph;
			last_h = h;
			
			//System.out.println(h+", "+ph);
			System.out.println(h+", "+ph);//+", "+delta_h_ph+", "+delta_p);
			
		}
		total_e /= Math.log(2);
		//System.out.println("E(H) = "+(total_h_ph/total_ph));
		System.out.println("E(H) = "+(total_h_ph));
		System.out.println("H(H) = "+(0-total_e));
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
		return y*y*(2.0*Math.log(y)-1.0) / (4.0*b);
	}
}
