package distributions.interfaces;


public interface PosteriorDistribution<THyper, TTheta, TData> {
	public void setPriorDistribution( PriorDistribution<THyper, TTheta, TData> prior);
	double[] getEntropyAndProbabilityAtTheta(TTheta thetas, TData[] data);
	double getProbabilityOfDataGivenTheta(TTheta thetas, TData[] data);
	public double getProbabilityOfThetaGivenData(TTheta thetas, TData[] data);
	public double getEntropyGivenTheta(TTheta theta);
	public double getDerivOfEntropyGivenTheta(TTheta thetas);
	public double[] getDerivsOfEntropyGivenTheta(TTheta thetas);
	public PriorDistribution<THyper, TTheta, TData> getConjugatePrior();
	public double getProbabilityOfDataGivenThetaPermuted(TTheta theta, TData[] data); 
}