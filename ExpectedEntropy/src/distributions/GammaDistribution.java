package distributions;  

import util.Functions;
import distributions.interfaces.PriorDistribution;

public class GammaDistribution implements PriorDistribution<double[],Double,Double> {
	double k;
	double o;
	double normalizing_constant = 1;

	public double getPriorProbabilityOfTheta(Double theta) {
		return normalizing_constant* Math.pow(theta,k-1)*Math.exp(-theta/o);
	}

	public void setHyperParameters(double[] hyperParameters) {
		k = hyperParameters[0];
		o = hyperParameters[1];
		normalizing_constant = 1.0/(Math.pow(o,k)*Functions.gamma(k));
	}
	
	public void updateFromData(Double[] data) {
		// TODO Auto-generated method stub
		
	}

}
