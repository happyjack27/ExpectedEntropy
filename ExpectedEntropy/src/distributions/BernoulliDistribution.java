package distributions;

import distributions.interfaces.*;

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
		return -(theta*Math.log(theta)+(1-theta)*Math.log(1-theta));
	}
	
	public double getDerivOfEntropyGivenTheta(Double theta) {
		return -Math.log(theta/(1-theta));
	}

	public double getProbabilityOfDataGivenTheta(Double theta, Integer[] data) {
		double zeros = data[0];
		double ones = data[1];
		return Math.pow(theta,zeros)*Math.pow(1-theta,ones);
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
