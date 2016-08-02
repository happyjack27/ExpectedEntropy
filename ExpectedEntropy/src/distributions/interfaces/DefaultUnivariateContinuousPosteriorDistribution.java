package distributions.interfaces;


public abstract class DefaultUnivariateContinuousPosteriorDistribution<TTheta> implements PosteriorDistribution<double[],TTheta,Double> {
	public boolean prior_is_on_H = false;
	public PriorDistribution<double[],TTheta,Double> prior;

	public DefaultUnivariateContinuousPosteriorDistribution() {
		super();
		setPriorDistribution(getConjugatePrior());
	}
	
	public DefaultUnivariateContinuousPosteriorDistribution(PriorDistribution<double[],TTheta,Double> prior) {
		super();
		setPriorDistribution(prior);
	}
	
	public void setPriorDistribution(PriorDistribution<double[],TTheta,Double> prior) {
		this.prior = prior;
	}
	
	public double[] getDerivsOfEntropyGivenTheta(TTheta thetas) {
		return new double[]{getDerivOfEntropyGivenTheta(thetas)};
	}

	
	public double[] getEntropyAndProbabilityAtTheta(TTheta thetas, Double[] data) {
		double h = getEntropyGivenTheta(thetas);
		double pDO = getProbabilityOfDataGivenThetaPermuted(thetas, data);
		double dHdO = getDerivOfEntropyGivenTheta(thetas);
		double pO = prior.getPriorProbabilityOfTheta(thetas);
		double ph = pDO*pO / ( prior_is_on_H ? 1 : Math.abs(dHdO));
		/*
		if( prior_is_on_H)
			pO = prior.getPriorProbabilityOfTheta((TTheta)h);
			*/
		
		return new double[]{h,ph,pDO,dHdO,pO};
	}

	public double getProbabilityOfThetaGivenData(TTheta theta, Double[] data) {
		return getProbabilityOfDataGivenTheta(theta,data)*prior.getPriorProbabilityOfTheta(/*prior_is_on_H ? getEntropyGivenTheta(theta) :*/ theta);
	}
}
