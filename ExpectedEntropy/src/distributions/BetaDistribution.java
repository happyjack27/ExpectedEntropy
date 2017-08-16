package distributions;

import distributions.interfaces.PriorDistribution;
import util.*;
import org.apache.commons.math3.util.*;

public class BetaDistribution implements PriorDistribution<double[],Double,Integer> {
	double[] hyperParameters;
	double multiplier = 1;
	double a = 1;
	double b = 1;
	
	BetaDistribution() {
		super();
		a = 1;
		b = 1;
		multiplier = Functions.gamma(a+b)/(Functions.gamma(a)*Functions.gamma(b));
	}
	BetaDistribution(double[] hyperParameters) {
		this();
		setHyperParameters(hyperParameters);
	}

	public double getPriorProbabilityOfTheta(Double theta) {
		return multiplier*FastMath.pow(theta,a-1)*FastMath.pow(1-theta,b-1);
	}

	public void setHyperParameters(double[] hyperParameters) {
		this.hyperParameters = hyperParameters;
		a = hyperParameters[0];
		b = hyperParameters[1];
		multiplier = Functions.gamma(a+b)/(Functions.gamma(a)*Functions.gamma(b));
	}
	
	public void updateFromData(Integer[] data) {
		a += data[0];
		b += data[1];
		multiplier = Functions.gamma(a+b)/(Functions.gamma(a)*Functions.gamma(b));
	}
}
