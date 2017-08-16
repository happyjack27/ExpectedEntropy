package distributions;

import distributions.interfaces.*;
import org.apache.commons.math3.util.*;


public class BernoulliDistribution extends DefaultUnivariateDiscretePosteriorDistribution<Double> implements hasMLE<Double,Integer> {
	
	public BernoulliDistribution() { super(); }
	public BernoulliDistribution(PriorDistribution<double[], Double, Integer> prior) { super(prior); }

	public PriorDistribution<double[],Double,Integer> getConjugatePrior() {
		return new BetaDistribution();
	}

	public double getEntropyGivenTheta(Double theta) {
		if( theta == 0 || theta == 1) {
			return 0;
		}
		return -(theta*FastMath.log(theta)+(1-theta)*FastMath.log(1-theta));
	}
	
	public double getDerivOfEntropyGivenTheta(Double theta) {
		return -FastMath.log(theta/(1-theta));
	}

	public double getProbabilityOfDataGivenTheta(Double theta, Integer[] data) {
		double zeros = data[0];
		double ones = data[1];
		return FastMath.pow(theta,zeros)*FastMath.pow(1-theta,ones);
	}
	
	public double getProbabilityOfDataGivenThetaPermuted(Double theta, Integer[] data) {
		return getProbabilityOfDataGivenTheta(theta,data)+getProbabilityOfDataGivenTheta(1-theta,data);
	}
	
	public Double getMLE(Integer[] data) {
		double q = 0;
		for( int i = 0; i < data.length; i++) {
			q += data[i];
		}
		return 1.0 - (q/(double)data.length);
	}
}
