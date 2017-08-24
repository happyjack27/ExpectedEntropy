package distributions;

import util.Functions;
import distributions.interfaces.PriorDistribution;
import org.apache.commons.math3.util.*;


public class DirichletDistribution implements PriorDistribution<double[],double[],Integer> {
	double[] a = new double[]{1,1};
	double normalizing_constant = 1;
	double log_normalizing_constant = 1;

	public double getLogPriorProbabilityOfTheta(double[] theta) {
		if( a.length != theta.length) {
			a = new double[theta.length];
			for( int i = 0; i < a.length; i++) {
				a[i] = 1;
			}
			setHyperParameters(a);
		}
		/*
		double m = 1;
		for( int i = 0; i < theta.length; i++) {
			m *= FastMath.pow(theta[i],a[i]-1);
		}
		return normalizing_constant*m;
		*/
		
		double m = 0;
		for( int i = 0; i < theta.length; i++) {
			if( theta[i] == 0) {
				continue;
				//return 0;
			}
			m += FastMath.log(theta[i])*(a[i]-1);
		}
		return log_normalizing_constant+m;
		//return normalizing_constant*FastMath.exp(m);
	}

	public double getPriorProbabilityOfTheta(double[] theta) {
		if( a.length != theta.length) {
			a = new double[theta.length];
			for( int i = 0; i < a.length; i++) {
				a[i] = 1;
			}
			setHyperParameters(a);
		}
		double m = 1;
		for( int i = 0; i < theta.length; i++) {
			m *= FastMath.pow(theta[i],a[i]-1);
		}
		return normalizing_constant*m;
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
		
		normalizing_constant =  Functions.gamma(bottom)/FastMath.exp(top);
		log_normalizing_constant =  FastMath.log(Functions.gamma(bottom)) - top;
	}
	

}
