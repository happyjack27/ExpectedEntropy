package distributions.interfaces;

import org.apache.commons.math3.util.*;

public abstract class DefaultUnivariateDiscretePosteriorDistribution<TTheta> implements PosteriorDistribution<double[],TTheta,Integer>  {

		public boolean prior_is_on_H = false;
		public PriorDistribution<double[],TTheta,Integer> prior;

		public DefaultUnivariateDiscretePosteriorDistribution() {
			super();
			setPriorDistribution(getConjugatePrior());
		}
		
		public DefaultUnivariateDiscretePosteriorDistribution(PriorDistribution<double[],TTheta,Integer> prior) {
			super();
			setPriorDistribution(prior);
		}
			
		public void setPriorDistribution(PriorDistribution<double[],TTheta,Integer> prior) {
			this.prior = prior;
		}
		
		public double[] getDerivsOfEntropyGivenTheta(TTheta thetas) {
			return new double[]{getDerivOfEntropyGivenTheta(thetas)};
		}

		public double[] getEntropyAndProbabilityAtTheta(TTheta thetas, Integer[] data) {
			double h = getEntropyGivenTheta(thetas);
			double pDO = getProbabilityOfDataGivenThetaPermuted(thetas, data);
			double dHdO = getDerivOfEntropyGivenTheta(thetas);
			double pO = prior.getPriorProbabilityOfTheta(thetas);
			double ph = pDO*pO / ( prior_is_on_H ? 1 : FastMath.abs(dHdO));
			//System.out.println("dHdO: "+dHdO+" pO: "+pO+" pDO: "+pDO+" ph: "+ph+" h: "+h+" p is on h: "+prior_is_on_H);
			/*
			if( prior_is_on_H)
				pO = prior.getPriorProbabilityOfTheta(h);
				*/
			
			return new double[]{h,ph,pDO,dHdO,pO};
		}

		public double getProbabilityOfThetaGivenData(TTheta theta, Integer[] data) {
			return getProbabilityOfDataGivenTheta(theta,data)*prior.getPriorProbabilityOfTheta(/*prior_is_on_H ? getEntropyGivenTheta(theta) :*/ theta);
		}
	}
