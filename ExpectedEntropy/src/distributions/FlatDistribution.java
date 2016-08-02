package distributions;

import distributions.interfaces.*;

public class FlatDistribution implements PriorDistribution<double[],Double,Double>, IntegrablePDF<Double,Double> {

	public double getCumulativeDistribution(Double thetas, Double value) {
		return value;
	}

	public double getPriorProbabilityOfTheta(Double theta) {
		return 1;
	}

	public void setHyperParameters(double[] hyperParameters) {
		
	}

	public void updateFromData(Double[] data) {
	}
}
