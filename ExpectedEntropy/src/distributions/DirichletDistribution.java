package distributions;

import util.Functions;
import distributions.interfaces.PriorDistribution;


public class DirichletDistribution implements PriorDistribution<double[],double[],Integer> {
	double[] a = new double[]{1,1};
	double normalizing_constant = 1;

	public double getPriorProbabilityOfTheta(double[] theta) {
		double m = 1;
		for( int i = 0; i < theta.length; i++) {
			m *= Math.pow(theta[i],a[i]-1);
		}
		return normalizing_constant*m;
		
		/*
		double m = 0;
		for( int i = 0; i < theta.length; i++) {
			if( theta[i] == 0) {
				continue;
				//return 0;
			}
			m += Math.log(theta[i])*(a[i]-1);
		}
		return normalizing_constant*Math.exp(m);
		*/
	}

	public void setHyperParameters(double[] hyperParameters) {
		a = hyperParameters;
		recalcNormalizingConstant();
		
	}

	public void updateFromData(Integer[] data) {
		for( int i = 0; i < data.length && i < a.length; i++) {
			a[i] += data[i];
		}
		recalcNormalizingConstant();
	}
	
	public void recalcNormalizingConstant() {
		double top = 0;
		double bottom = 0;
		for( int i = 0; i < a.length; i++) {
			top += Functions.logGamma(a[i]);
			bottom += a[i];
		}
		
		normalizing_constant =  Functions.gamma(bottom)/Math.exp(top);
	}
	

}
