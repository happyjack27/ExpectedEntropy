package distributions.interfaces;

public interface PriorDistribution<THyper, TTheta, TData> { //todo, generalize by adding <data,param> to definition
	public double getPriorProbabilityOfTheta(TTheta theta);
	public void setHyperParameters(THyper hyperParameters);
	public void updateFromData(TData[] data);
}