package distributions;

import java.util.*;

import distributions.interfaces.PosteriorDistribution;
import distributions.interfaces.PriorDistribution;
import util.*;
import org.apache.commons.math3.*;
import org.apache.commons.math3.distribution.*;

// http://www.sortie-nd.org/lme/Statistical%20Papers/Burnham_and_Anderson_2004_Multimodel_Inference.pdf
// http://tuvalu.santafe.edu/~simon/page7/page9/page9.html

public class CategoricalDistribution implements PosteriorDistribution<double[],double[],Integer> {

	public boolean prior_is_on_H = false;
	public boolean use_prior_for_sampling = false;

	double[] thetas;
	int[][] permutations;
	
	PriorDistribution<double[],double[],Integer> prior;

	public static int add_to_all_cats = 0;
	
	public static double GAMMA = 0.577215664901532860606512090082;
	public static double GAMMA_MINX = 1.e-12;
	public static double DIGAMMA_MINNEGX = -1250;
	public static double C_LIMIT = 49;
	public static double S_LIMIT = 1e-5;
	
	
	/*
	
		BEST PAPERS EVER!!!
	
	//  http://dukespace.lib.duke.edu/dspace/bitstream/handle/10161/2458/D_Paisley_John_a_201005.pdf?sequence=1
	//  http://www.jmlr.org/papers/volume15/archer14a/archer14a.pdf
	 * 
	 * 
	 * 
	 */
	
	
	/*
	 * 
	 * 
	 * --- maybe just apply the loss function to h right away, and integrate there.
	 * for instance, lets say the loss function is just h, then just take integral( h*dx) / integral(dx)
C=e^{H=-\sum_{i \in I}p_iln(p_i)} = \prod_{i \in I}p_i^{-p_i},
p(P|A) \propto_{f(A))} \prod_{i \in I}p_i^{a_i},

E(C|A) =
\iint_{\sum_{i \in I}p_i = 1} C \cdot p(C|A) \cdot dp
\propto_{f(A)}
\iint_{\sum_{i \in I}p_i = 1} 
\prod_{i \in I} p_i^{a_i-p_i}
\cdot dp_i

C=e^{H=-\sum_{i \in I}p_iln(p_i)} = \prod_{i \in I}p_i^{-p_i},
p(P|A) \propto \prod_{i \in I}p_i^{a_i},

E(C|A) \propto_{a}
\iint_{\sum_{i \in I}p_i = 1} C \cdot p(C|A) \cdot dp
=
\iint_{\sum_{i \in I}p_i = 1} 
\prod_{i \in I} p_i^{a_i-p_i}
\cdot dp_i


======

E(H) = 

E(H) = 
\iint_{\sum_{i \in I}p_i = 1} 
(\sum_{i \in I} p_i ln(p_i))
\cdot \prod_{i \in I} p_i^{a_i} 
\cdot dp_i

-or-

E(H) = 
\iint_{\sum_{i \in I}p_i = 1} 
(\sum_{i \in I} p_i ln(p_i))
\cdot \prod_{i \in I} p_i^{a_i} 
\cdot (ln(p_i)+1) 
\cdot dp_i

integral{ 
-sum_I (pi ln pi)
*product_I(pi^ai)
*product_I(dpi)
?/?*(|product_I(dH/dp)=(lnpi+1)|)
} over the region
p0 = 0 ... 1
p1 = 0 ... 1-p0
p2 = 0 ... 1-p0-p1
p3 = 0 ... 1-p0-p1-p2
pn = 1-p0-p1-p2-...p(n-1)

/ Dirichlet Beta(ai's..)


C=e^{H=\sum_{i \in I}p_iln(p_i)}, 
E(C) = 
\iint_{\sum_{i \in I}p_i = 1} 
\prod_{i \in I} p_i^{a_i-p_i}
\cdot dp_i
- or -

E(C) = 
\iint_{\sum_{i \in I}p_i = 1} 
\prod_{i \in I} p_i^{a_i-2p_i}
\cdot (ln(p_i)+1)
\cdot dp_i

integral{ 
product_I (e^(-pi ln pi))
*product_I(pi^ai)
*product_I(dpi)
?/?*(|product_I(dC/dp)=(-(ln(pi)+1)*e^(-pi-ln(pi)))|)
} over the region
p0 = 0 ... 1
p1 = 0 ... 1-p0
p2 = 0 ... 1-p0-p1
p3 = 0 ... 1-p0-p1-p2
pn = 1-p0-p1-p2-...p(n-1)

/ Dirichlet Beta(ai's..)



p(X=x) = x^(counts(a)*(1-x)^(counts(not a)) / constant func of counts (B(a+1,not a+1))
 = e^(lnx*counts(a))* e^(ln(1-x)*counts(not a))
===========
inverse of I = y=-xlnx , x= -y/W(-y) or y=-lnx, e^-y=x

----------
wright function:

x = yln(y) 

y = x * w(lnx) = x/W(x) = e^W(x)

dx/dy = lny + 1

dy/dx = 1/(lny+1)  = 1/(ln(x*w(lnx))+1) = 1/(ln(x)+ln(w(lnx))+1) = 1/(W(x)+1)

----------

y = number of categories the selection represents.
x = probability of the selection.

inverse of # categories = y = e^(-xlnx) , x = -lny/W(-lny) or y=1/x, x=1/y

dy/dx = -(ln(x)+1)*e^(-x-ln(x)) = -(W(xlnx)+1)*e^(-xlnx) = -y*(W(-log(y))+1)

dx/dy = -1 / [ y*(W(-log(y))+1) ]

--------------

= -[ lny * w'(-log(y)) + W(-log(y)) ] / [y * W(-log(y))^2 ]

= -[ lny * [ W(-log(y))/(-log(y)[1+W(-log(y))]) ] + W(-log(y)) ] / [y * W(-log(y))^2 ]


...so... does (-ln(x)-1)e^(-x*ln(x))* = y*(W(-log(y))+1) + C?  when y=e^(-xlnx)

y*(W(-log(y))+1) = e^(-xlnx)*(W(-log(e^(-xlnx)))+1) = e^(-xlnx)*(W(xlnx)+1)

so does W(xlnx)+1 = -ln(x)-1?

ln(x)*e^(lnx)=xlnx


W'(x) = W(x)/(x[1+W(x)])

[ W(-log(y))/(-log(y)[1+W(-log(y))]) ]

=============

	 */
	
	
	//=================
	
	
	/*
number of coins remainder must be = (target h - current h)/(1-current x) = (H-h)/(1-x)

change in # coins needed per change in p used:

d'/dx = (H-h) / (1-x)^2

change in # coins needed per change in h created:

d'/dh = - 1 / (1-p)


number of categories needed = e^((H-h)/(1-x))

d'/dx = (H-h)e^[(H-h)/(1-x)] / (1-x)^2

d'/dh = - e^[(H-h)/(1-p)] / (1-p)

number of categories used = e^H.

is there some kind of recursion relationship?

number of ways to pick p cats from n is n:p = gama(n+1)/[gama(p+1)gamma(p-n+1)]

# categories a thing represents = 1/p.
... / (1-p)


--

one coin, needs 1 bit:
first category:
1/2 a bit is done, with 1/2 the probability used up.
the remaining category must pick up the remaining 1/2 a bit.
conditioned on the previous categories not being chosen, so p(a|~b) = p(a)*p(~b|a)/p(~b) = 1 (when down to 2) e^(1ln2/(0.5)) = e^2ln2 = 2

the probability of the last bit must be 1/(number of categories needed after the bit before)

--

3 categories, needs x bits:
first category:
h a bit is done, with p the probability used up.
the remaining 2 categories must pick up the remaining x-h bits.

to do so they need e^((x-h)ln2/(1-p)) categories among them, given that it is known that one of them is selected.

so each probability is divided by (1-p).

e^[[-p(x4)*ln(p(x4))-p(x5)*ln(p(x5))...]/p(!x1&!x2&!x3..)] = target categories - number covered earlier.

now we have a weighted coin flip, where the weights dont need to add up to 1 (they add up to 1-p).

last category must have probability 1-p-q, which must equal 1/(number of categories needed after the bit before)

  1-p-q = 1/(e^((x-h)ln2/(1-p)) - q/(1-p) * ln(q))

..solve for q


p(b|~a)*1/info = target - p(a)*I(a)


integral domain, where X = number of categories needed, N = total # categories, n = current category ( left to right)

 0 - lessor of( N-n,...)
 
 ==== completely different way:
 
 for each unused category, calculate the likilihood that it will not be used.  since its either used or not. and we only care if its certainly never used.
 if we're wrong - the wrongness can be justified by saying "there's some other latent state we didn't know about, that flipped." 

=======

equal H is reached when 
s (change in H must equal 0 ) and sum ( dO ) = 0 (sum of thetas must always = 1)

so solve... find vector x such that sum(x*y)=0=sum(x) where y=dHdO vector.  since we're keeping sum dO=0, we can use just the d(plnp)=lnp+1 part? 

sum(y) = 0;
sum(ylnx+y)=0;

================

	 */
	
	
	/*
	 * 
	 * dx - distribute counterbalance on other probabilities, as if adding a count to one. - subtract it to other in proportion to their probability
	 * 
	 * [(a+h) ln (a+h) + (b-h*b/(1-a)) * ln (b-h*b/(1-a)) + (c-h*c/(1-a)) * ln (c-h*c/(1-a)) - (alna+blnb+clnc)] / h...  lim h > 0 
	 * 
	 * = [(1-a+h)*ln(1-a+h) + b*(1-h/a)*ln(b*(1-h/a)) + c*(1-h/a) * ln(c*(1-h/a)) 
	 * - ((1-a)ln(1-a)+blnb+clnc)] / h.
	 * 
	 * = [(1-a+h)*ln(1-a+h) + b*(1-h/a)*[ln(b) + ln(1-h/a)] + c*(1-h/a) * (ln(c)+ln(1-h/a)) 
	 * - ((1-a)ln(1-a)+blnb+clnc)] / h.
	 * 
	 * = [(1-a+h)*ln(1-a+h) + (1-h/a) * [b*ln(b) + c*ln(c)] + (b+c=a)*(1-h/a) * ln(1-h/a)
	 * - ((1-a)ln(1-a)+blnb+clnc)] / h.
	 * 
	 * = [(1-a+h)*ln(1-a+h) + (1-h/a-1) * [b*ln(b) + c*ln(c)] + (a-h)*ln(1-h/a) - (1-a)ln(1-a)] / h.
	 * 
	 * = [(1-a+h)*ln(1-a+h) + (a-h)*ln(1-h/a) - (1-a)ln(1-a) - (h/a) * [b*ln(b) + c*ln(c)]] / h.
	 * 
	 * = [(1-a)*ln(1-a+h) + h*ln(1-a+h) + a*ln(1-h/a) - h*ln(1-h/a) - (1-a)ln(1-a) - (h/a) * [b*ln(b) + c*ln(c)]] / h.

	 * = [(1-a)*ln(1-a+h) + a*ln(1-h/a) - (1-a)ln(1-a)] / h  + ln(1-a+h) - ln(1-h/a) - (1/a) * [b*ln(b) + c*ln(c)]

	 * = a*(ln(a-h)-ln(a)) / h + (a-1)(ln(1-a)-ln(1-a+h))/h  + ln(1-a+h) - ln(1-h/a) - (1/a) * [b*ln(b) + c*ln(c)]
	 * 
	 * = -a*(1/a) + (1-a)*(1/(1-a)) + ...
	 * 
	 * = -1  + 1 + ln(1-a+h) - ln(1-h/a) - (1/a) * [b*ln(b) + c*ln(c)]
	 * 
	 * = ln(1-a) - ln(1) - (1/a) * [b*ln(b) + c*ln(c)]
	 * 
	 * d(entropy)/d(1-a) = ln(1-a) - [b*ln(b) + c*ln(c)]/a
	 * 
	 * where a = b+c+....
	 * 
	 * try just estimating the information for that one choice, given the probability of the choice, given it's chosen:  
	 * that is, estimate the density function for log(p)=x  inverse function is of course e^x.
	 *
	 * 
	 * -- then we can do conditional distribution for moving to the next (that is, given that the likilihood of symbol a is x)
	 */

	public void setNumberOfCategories(int n) {
		thetas = new double[n];
		double[] hp = new double[n];
		for( int i = 0; i < hp.length; i++) {
			hp[i] = 1;
		}
		prior.setHyperParameters(hp);
		
	}

	public CategoricalDistribution() {
		super();
		setPriorDistribution(getConjugatePrior());
	}
	
	public CategoricalDistribution(PriorDistribution<double[],double[],Integer> prior) {
		super();
		setPriorDistribution(prior);
	}
		
	
	public void setPriorDistribution(PriorDistribution<double[], double[], Integer> prior) {
		this.prior = prior;
		
	}
	public Integer[] jointToMarginData(Integer[] data,int[][] marginal_categories) {
		Integer[] ii = new Integer[marginal_categories.length];
		for( int i = 0; i < ii.length; i++) {
			ii[i] = 0;
			for( int k : marginal_categories[i]) {
				ii[i] += data[k];
			}
		}
		return ii;
	}
	public double[] jointToMarginTheta(double[] data,int[][] marginal_categories) {
		double[] ii = new double[marginal_categories.length];
		for( int i = 0; i < ii.length; i++) {
			ii[i] = 0;
			for( int k : marginal_categories[i]) {
				ii[i] += data[k];
			}
		}
		return ii;
	}
	
	//if you want the conditional entropy of X given Y, the marginal categories should be the categories for Y.
	public double getConditionalEntropy(Integer[] data,int[][] marginal_categories, int num_samples) {
		Integer[] margin_data = jointToMarginData(data,marginal_categories);
		if( data.length == 0) {
			return 0;
		}
		
		//NEED TO GO THROUGH THIS
		setNumberOfCategories(data.length);

		Vector<double[]> results = new Vector<double[]>();
		Vector<Pair<Double,double[]>> samples = new Vector<Pair<Double,double[]>>();
		
		//do a bunch of random samples
		for( int i = 0; i < num_samples; i++) {		
			//calculate entropy and probability
			double[] thetas = getRandomThetas(data);
			double[] margin_thetas = jointToMarginTheta(thetas,marginal_categories);
			
			double[] dd = getEntropyAndLogProbabilityAtTheta(thetas,data,false);
			double[] dd_margin = getEntropyAndLogProbabilityAtTheta(margin_thetas,margin_data,false);
			dd[0] -= dd_margin[0];
			
			samples.add(new Pair<Double,double[]>(dd[0],dd));
		}
		double[] thetas = getMLEThetas(data);
		double[] thetas2 = getMLEThetas(margin_data);
		double[] dd2 = getEntropyAndLogProbabilityAtTheta(thetas,data,false);
		double[] dd2_m = getEntropyAndLogProbabilityAtTheta(thetas2,margin_data,false);
		double[] dd0 = getEntropyAndLogProbabilityAtTheta(getMinEntropyThetas(data.length),data,false);
		double[] dd0_m = getEntropyAndLogProbabilityAtTheta(getMinEntropyThetas(margin_data.length),margin_data,false);
		double[] dd1 = getEntropyAndLogProbabilityAtTheta(getMaxEntropyThetas(data.length),data,false);
		double[] dd1_m = getEntropyAndLogProbabilityAtTheta(getMaxEntropyThetas(margin_data.length),margin_data,false);
		dd0[0] -= dd0_m[0];
		dd1[0] -= dd1_m[0];
		dd2[0] -= dd2_m[0];
		samples.add(new Pair<Double,double[]>(dd0[0],dd0));
		samples.add(new Pair<Double,double[]>(dd1[0],dd1));
		samples.add(new Pair<Double,double[]>(dd2[0],dd2));
		//getEntropyAndProbabilityAtTheta(getMinEntropyThetas(),data,false);
		//doulbe[] min = 
		Collections.sort(samples);
		//results.add
		
		//adjust log p arithmetically to max of 64, then take the exponent;
		double max_log_p = samples.get(0).b[1];
		for(int i = 0; i < samples.size(); i++) {
			if( samples.get(i).b[1] != samples.get(i).b[1]) {
				continue;
			}
			if( samples.get(i).b[1] > max_log_p || max_log_p != max_log_p) {
				max_log_p = samples.get(i).b[1];
			}
		}
		for(int i = 0; i < samples.size(); i++) {
			if( samples.get(i).b[1] != samples.get(i).b[1]) {
				samples.get(i).b[1] = 0;
			}
			samples.get(i).b[1] = Math.exp(samples.get(i).b[1] + 64.0 - max_log_p);
		}
		
		//now add all samples to results.
		for(int i = 0; i < samples.size(); i++) {
			results.add(samples.get(i).b);
		}
		double maxH = Math.log((double)data.length);/////Math.log(2);
		
		//integrate
		double sumP = 0;
		double sumHP = 0;
		double last_h = results.get(0)[0];
		double last_p = results.get(0)[1];

		for( int i = 1; i < results.size(); i++) {
			double h = results.get(i)[0];
			double p = results.get(i)[1];
			double dH = h-last_h;
			double avgH = (h+last_h)/2.0;
			double avgP = (p+last_p)/2.0;
			double dHP = avgH*dH*avgP;
			double dP = dH*avgP;
			if( dP == 0) {
				//System.out.println("dP is 0 "+dH+" "+avgP);
			}
			if( dP == dP && dHP == dHP) {
				sumHP += dHP;
				sumP += dP;
			} else {
				//System.out.println("NAN! :"+dP+" "+dHP+" | "+avgH+" "+dH+" "+avgP);
			}
			last_h = h;
			last_p = p;
		}
		//System.out.println("sums "+sumHP+" "+sumP);

		return (sumHP/sumP);/////Math.log(2);
	}
	public double getAtPercentile(Vector<double[]> results, double percentile) {
		if( results.size() == 0) {
			return 0;
		}
		return getAtSortedPercentiles(results, new double[]{percentile})[0];
	}
	public double[] getAtSortedPercentiles(Vector<double[]> results, double[] percentiles) {
		double[] ret = new double[percentiles.length];
		//integrate
		double sumP = 0;
		double last_h = results.get(0)[0];
		double last_p = results.get(0)[1];
		for( int i = 1; i < results.size(); i++) {
			double h = results.get(i)[0];
			double p = results.get(i)[1];
			double dH = h-last_h;
			double avgP = (p+last_p)/2.0;
			double dP = dH*avgP;
			if( dP == dP) {
				sumP += dP;
			} else {
				//System.out.println("NAN! :"+dP+" "+dHP+" | "+avgH+" "+dH+" "+avgP);
			}
			last_h = h;
			last_p = p;
		}
		double totP = sumP;
		double targetP = percentiles[0] * totP;

		int n = 0;
		sumP = 0;
		last_h = results.get(0)[0];
		last_p = results.get(0)[1];
		for( int i = 1; i < results.size(); i++) {
			double h = results.get(i)[0];
			double p = results.get(i)[1];
			double dH = h-last_h;
			double avgP = (p+last_p)/2.0;
			double dP = dH*avgP;
			if( dP == dP) {
				sumP += dP;
			} else {
				//System.out.println("NAN! :"+dP+" "+dHP+" | "+avgH+" "+dH+" "+avgP);
			}
			while( sumP > targetP) {
				double prevP = sumP - dP;
				ret[n] = (last_h + dH*(targetP - prevP)/dP);/////Math.log(2);
				n++;
				if( n >= percentiles.length) {
					break;
				}
				targetP = percentiles[n] * totP;
				//return (h+last_h)/2.0;
			}
			if( n >= percentiles.length) {
				break;
			}
			last_h = h;
			last_p = p;
		}
		
		return ret;
	}
	
	public double integrateDiscrete(Vector<double[]> results) {
		if( results.size() == 0) {
			return 0;
		}
		
		//integrate
		double sumP = 0;
		double sumHP = 0;
		for( int i = 0; i < results.size(); i++) {
			double h = results.get(i)[0];
			double p = results.get(i)[1];
			sumHP += h*p;
			sumP  += p;
		}
		return (sumHP/sumP);/////Math.log(2);
	}	
	public double integrateSortedEntropyCurve(Vector<double[]> results) {
		if( results.size() == 0) {
			return 0;
		}
		
		//integrate
		double sumP = 0;
		double sumHP = 0;
		double last_h = results.get(0)[0];
		double last_p = results.get(0)[1];

		for( int i = 1; i < results.size(); i++) {
			double h = results.get(i)[0];
			double p = results.get(i)[1];
			double dH = h-last_h;
			double avgH = (h+last_h)/2.0;
			double avgP = (p+last_p)/2.0;
			double dHP = avgH*dH*avgP;
			double dP = dH*avgP;
			if( dP == 0) {
				//System.out.println("dP is 0 "+dH+" "+avgP);
			}
			if( dP == dP && dHP == dHP) {
				sumHP += dHP;
				sumP += dP;
			} else {
				//System.out.println("NAN! :"+dP+" "+dHP+" | "+avgH+" "+dH+" "+avgP);
			}
			last_h = h;
			last_p = p;
		}
		//System.out.println("sums "+sumHP+" "+sumP);

		return (sumHP/sumP);/////Math.log(2);
	}	
	public double getEntropyOfSortedEntropyCurve(Vector<double[]> results) {
		double res = 10.0/(double)results.size();//Math.exp(-Math.log((double)results.size())*0.5);
		double maxH = results.get(results.size()-1)[0];///Math.log(2); 
		//System.out.println("maxH: "+maxH);
		
		double[] percentiles = new double[(int)Math.ceil((double)1.0 / res)];
		for( int i = 0; i < percentiles.length; i++) {
			percentiles[i] = res*(double)i;
			if( percentiles[i] > 1) {
				percentiles[i] = 1;
			}
		}
		double[] hs = this.getAtSortedPercentiles(results, percentiles);
		
		double e = 0;
		double max_e = 0;
		for( int i = 1; i < hs.length; i++) {
			max_e -= res * Math.log(res);
		}
		for( int i = 1; i < hs.length; i++) {
			double dH = (hs[i] - hs[i-1])/maxH;
			e -= res * Math.log(dH);
		}
		return (e-max_e);/////Math.log(2);
	}	
	
	//if you want the conditional entropy of X given Y, the marginal categories should be the categories for Y.
	//H(X,Y),H(X),H(Y),H(X|Y),H(Y|X),MI(X,Y)
	public Vector<double[]>[] getAllEntropiesVectors(Integer[] data, int[][] marginal_categories1, int[][] marginal_categories2, int num_samples) {
		Vector<double[]>[]  rets = new Vector[6];
		if( data.length == 0) {
			return rets;
		}
		Integer[] margin_data1 = jointToMarginData(data,marginal_categories1);
		Integer[] margin_data2 = jointToMarginData(data,marginal_categories2);
	
		setNumberOfCategories(data.length);
	
		Vector<Pair<Double,double[]>> vhxy = new Vector<Pair<Double,double[]>>();
		Vector<Pair<Double,double[]>> vhx = new Vector<Pair<Double,double[]>>();
		Vector<Pair<Double,double[]>> vhy = new Vector<Pair<Double,double[]>>();
		Vector<Pair<Double,double[]>> vhxgy = new Vector<Pair<Double,double[]>>();
		Vector<Pair<Double,double[]>> vhygx = new Vector<Pair<Double,double[]>>();
		Vector<Pair<Double,double[]>> vmi = new Vector<Pair<Double,double[]>>();
	
		Vector<Vector<Pair<Double,double[]>>> vs = new Vector<Vector<Pair<Double,double[]>>>();
		vs.add(vhxy);
		vs.add(vhx);
		vs.add(vhy);
		vs.add(vhxgy);
		vs.add(vhygx);
		vs.add(vmi);
	
		//do a bunch of random samples
		for( int i = 0; i < num_samples; i++) {
		//calculate entropy and probability
			double[] thetas = getRandomThetas(data);
			switch(i) {
				case 0:
				thetas = getMinEntropyThetas(data.length);
			break;
				case 1:
				thetas = getMaxEntropyThetas(data.length);
			break;
				case 2:
				thetas = getMLEThetas(data);
			break;
			}
			double[] margin_thetas1 = jointToMarginTheta(thetas,marginal_categories1);
			double[] margin_thetas2 = jointToMarginTheta(thetas,marginal_categories2);
		
			double[] hxy = getEntropyAndLogProbabilityAtTheta(thetas,data,false);
			double[] hx = getEntropyAndLogProbabilityAtTheta(margin_thetas1,margin_data1,false);
			double[] hy = getEntropyAndLogProbabilityAtTheta(margin_thetas2,margin_data2,false);
			double[] hxgy = hxy.clone();
			double[] hygx = hxy.clone();
			double[] mi = hxy.clone();
		
			hxgy[0] -= hy[0];
			hygx[0] -= hx[0];
			mi[0] = hx[0]+hy[0]-hxy[0];
		
			vhxy.add(new Pair<Double,double[]>(hxy[0],hxy));
			vhx.add(new Pair<Double,double[]>(hx[0],hx));
			vhy.add(new Pair<Double,double[]>(hy[0],hy));
			vhxgy.add(new Pair<Double,double[]>(hxgy[0],hxgy));
			vhygx.add(new Pair<Double,double[]>(hygx[0],hygx));
			vmi.add(new Pair<Double,double[]>(mi[0],mi));
		}
	
		//double[] rets = new double[vs.size()];
		for( int vi = 0; vi < rets.length; vi++) {
			Vector<double[]> results = new Vector<double[]>();
			Vector<Pair<Double,double[]>> vp = vs.get(vi);
			Collections.sort(vp);
			//adjust log p arithmetically to max of 0, then take the exponent;
			double max_log_p = vp.get(0).b[1];
			for(int i = 0; i < vp.size(); i++) {
				if( vp.get(i).b[1] > max_log_p) {
					max_log_p = vp.get(i).b[1];
				}
			}
		
			for(int i = 0; i < vp.size(); i++) {
				vp.get(i).b[1] = Math.exp(vp.get(i).b[1] + 0.0 - max_log_p);
			}
		
			//now add all samples to results.
			for(int i = 0; i < vp.size(); i++) {
				results.add(new double[]{vp.get(i).b[0],vp.get(i).b[1]});
			}
		
			rets[vi] = results;
		}
	
		return rets;
	}

	public double getBayesianActualEntropy(Integer[] ii, double[] actual_dist, int num_samples) {
		if( num_samples == 0) { num_samples = 100; }
		
		Vector<double[]> results = new Vector<double[]>();
		Vector<Pair<Double,double[]>> samples = new Vector<Pair<Double,double[]>>(); 
		
		for( int i = 0; i < num_samples; i++) {
			double[] thetas = getRandomThetas(ii.length,false);
			if( i == 0) {
				//thetas = getMaxEntropyThetas(ii.length);
			}
			if( i == 1) {
				thetas = actual_dist.clone();
			}
			double logp_theta_given_data = 0;
			boolean ok = true;
			for( int j = 0; j < thetas.length; j++) {
				if( thetas[j] != 0) {
					logp_theta_given_data += ((double)ii[j]) * Math.log(thetas[j]);
				} else if( ii[j] != 0) {
					ok = false;
					break;
				}
			}
			if( !ok) { //p=0
				continue;
			}
			
			logp_theta_given_data = Math.exp(logp_theta_given_data);
			
			for( int j = 0; j < thetas.length; j++) {
				if( thetas[j] != 0) {
					double H = -Math.log(thetas[j]);
					double p = actual_dist[j] * logp_theta_given_data;
					samples.add(new Pair<Double,double[]>(H,new double[]{H,p}));
				}
			}
		}
		Collections.sort(samples);
		//results.add
		for(int i = 0; i < samples.size(); i++) {
			results.add(samples.get(i).b);
		}		
		
		return integrateDiscrete(results);
		//return integrateSortedEntropyCurve(results);
	}
	
	public double getBayesianActualEntropy_old(Integer[] ii, double[] actual_dist, int samples) {
		if( samples == 0) { samples = 100; }
		
		Vector<double[]> hpVect = new Vector<double[]>();
		
		double total_kldp = 0;
		double total_p = 0;
		for( int i = 0; i < samples; i++) {
			double[] thetas = getRandomThetas(ii.length,false);
			double kld = 0;
			double lp = 0;
			for( int j = 0; j < thetas.length; j++) {
				if( thetas[j] != 0) {
					kld += -actual_dist[j] * Math.log(thetas[j]);//(Math.log(thetas[j]) - Math.log(actual_dist[j])); //can't do kldiv because using different model would be lower entropy
					lp += ((double)ii[j]) * Math.log(thetas[j]);
				} else {
					if( ii[j] > 0) {
						// then =*0^(x<>0), so p = 0.
						//otherwise, it's =*1.
						continue;
					} else {
						//0^a*p*0^0 / p*0^0 
						//what to do?
						
						//numerator is : numerator + a * log(0) * 0^0
						//denominator = 0
					}
				}
			}
			double p = Math.exp(lp);
			total_kldp += p*kld;
			total_p += p;
		}
		
		return total_kldp/total_p;
	}
	
	public double[] getMLEThetas(Integer[] data) {
		double[] thetas = new double[data.length];
		double tot = 0;
		for( int i = 0; i < data.length; i++) {
			tot += data[i];
		}
		if( tot == 0) {
			double p = 1.0/(double)data.length;
			for( int i = 0; i < data.length; i++) {
				thetas[i] = p;
			}
		} else {
			double p = 1.0/tot;
			for( int i = 0; i < data.length; i++) {
				thetas[i] = p * (double) data[i];
			}
		}
		// TODO Auto-generated method stub
		return thetas;
	}
	public Vector<double[]> getEntropyCurveLogarithmic(Integer[] data, int num_samples) {
		if( data.length == 0) {
			return new Vector<double[]>();
		}
		setNumberOfCategories(data.length);
	
		Vector<Pair<Double,double[]>> vhxy = new Vector<Pair<Double,double[]>>();
	
		//do a bunch of random samples
		for( int i = 0; i < num_samples; i++) {
		//calculate entropy and probability
			double[] thetas = getRandomThetas(data);
			switch(i) {
				case 0:
				thetas = getMinEntropyThetas(data.length);
			break;
				case 1:
				thetas = getMaxEntropyThetas(data.length);
			break;
				case 2:
				thetas = getMLEThetas(data);
			break;
			}
			double[] hxy = getEntropyAndLogProbabilityAtTheta(thetas,data,false);
		
		
			vhxy.add(new Pair<Double,double[]>(hxy[0],hxy));
		}
	
		Vector<double[]> results = new Vector<double[]>();
		Vector<Pair<Double,double[]>> vp = vhxy;
		Collections.sort(vp);
		//adjust log p arithmetically to max of 0, then take the exponent;
		double max_log_p = vp.get(0).b[1];
		for(int i = 0; i < vp.size(); i++) {
			if( vp.get(i).b[1] > max_log_p) {
				max_log_p = vp.get(i).b[1];
			}
		}
	
		for(int i = 0; i < vp.size(); i++) {
			vp.get(i).b[1] = Math.exp(vp.get(i).b[1] + 0.0 - max_log_p);
		}
	
		//now add all samples to results.
		for(int i = 0; i < vp.size(); i++) {
			results.add(new double[]{vp.get(i).b[0],vp.get(i).b[1]});
		}
		
		return results;
	}
	
	
	
	public Vector<double[]> getEntropyCurveMultiplicative(Integer[] data, int num_samples) {
		if( data.length == 0) {
			return new Vector<double[]>();
		}

		setNumberOfCategories(data.length);
		Integer[] data2;
		if( add_to_all_cats == 0) {
			data2 = data;
		} else {
			data2 = new Integer[data.length];
			for( int i = 0; i < data.length; i++) {
				data2[i] = data[i] + add_to_all_cats;
			}
			
		}
		
		Vector<double[]> results = new Vector<double[]>();
		Vector<Pair<Double,double[]>> samples = new Vector<Pair<Double,double[]>>();
		
		//do a bunch of random samples
		for( int i = 0; i < num_samples; i++) {		
			//calculate entropy and probability
			double[] dd = getEntropyAndProbabilityAtTheta(getRandomThetas(data2),data2,false);
			if( i == 0) {
				dd = getEntropyAndProbabilityAtTheta(getMinEntropyThetas(data2.length),data2,false);
			}
			if( i == 1) {
				dd = getEntropyAndProbabilityAtTheta(getMaxEntropyThetas(data2.length),data2,false);
			}
			if( i == 2) {
				dd = getEntropyAndProbabilityAtTheta(getMLEThetas(data2),data2,false);
			}
			samples.add(new Pair<Double,double[]>(dd[0],dd));
		}
		//getEntropyAndProbabilityAtTheta(getMinEntropyThetas(),data,false);
		//doulbe[] min = 
		Collections.sort(samples);
		//results.add
		for(int i = 0; i < samples.size(); i++) {
			results.add(samples.get(i).b);
		}
		return results;
	}

	public double[] getSummaryStats(Integer[] data, int num_samples) {
		if( data.length == 0) {
			return new double[]{0};
		}

		setNumberOfCategories(data.length);
		Integer[] data2;
		if( add_to_all_cats == 0) {
			data2 = data;
		} else {
			data2 = new Integer[data.length];
			for( int i = 0; i < data.length; i++) {
				data2[i] = data[i] + add_to_all_cats;
			}
			
		}
		
		Vector<double[]> results = new Vector<double[]>();
		Vector<Pair<Double,double[]>> samples = new Vector<Pair<Double,double[]>>();
		
		//do a bunch of random samples
		for( int i = 0; i < num_samples; i++) {		
			//calculate entropy and probability
			double[] dd = getEntropyAndProbabilityAtTheta(getRandomThetas(data2),data2,false);
			if( i == 0) {
				dd = getEntropyAndProbabilityAtTheta(getMinEntropyThetas(data2.length),data2,false);
			}
			if( i == 1) {
				dd = getEntropyAndProbabilityAtTheta(getMaxEntropyThetas(data2.length),data2,false);
			}
			if( i == 2) {
				dd = getEntropyAndProbabilityAtTheta(getMLEThetas(data2),data2,false);
			}
			samples.add(new Pair<Double,double[]>(dd[0],dd));
		}
		//getEntropyAndProbabilityAtTheta(getMinEntropyThetas(),data,false);
		//doulbe[] min = 
		Collections.sort(samples);
		//results.add
		for(int i = 0; i < samples.size(); i++) {
			results.add(samples.get(i).b);
		}
		double maxH = Math.log((double)data.length);/////Math.log(2);
		
		//integrate
		double sumP = 0;
		double sumHP = 0;
		double last_h = results.get(0)[0];
		double last_p = results.get(0)[1];

		for( int i = 1; i < results.size(); i++) {
			double h = results.get(i)[0];
			double p = results.get(i)[1];
			double dH = h-last_h;
			double avgH = (h+last_h)/2.0;
			double avgP = (p+last_p)/2.0;
			sumHP += avgH*dH*avgP;
			sumP += dH*avgP;
			last_h = h;
			last_p = p;
		}
		return new double[]{(sumHP/sumP)/Math.log(2.0)};
	}
	
	public static double[] getSummaryStatsForEntropyPDF( Vector<double[]> samples) {

		double last_h = 0;
		double last_ph = 0;
		double total_h_ph = 0;
		double total_ph = 0;
		double total_e = 0;
		
		//compute normalizing constant
		for( int i = 0; i < samples.size(); i++) {
			double[] result = samples.get(i);
			double h = result[0];///Math.log(2);
			double ph = result[1];
			double delta_h = h-last_h;
			double delta_ph = ph-last_ph;
			double delta_h_ph = 0;
			double delta_p = 0;
			delta_p = delta_h*(last_ph+ph)/2.0;
			total_ph += delta_p;
			last_ph = ph;
			last_h = h;
		}
		
		last_h = 0;
		last_ph = 0;
		for( int i = 0; i < samples.size(); i++) {
			double[] result = samples.get(i);
			double h = result[0];///Math.log(2);
			double ph = result[1]/total_ph;
			double delta_h = h-last_h;
			double delta_ph = ph-last_ph;
			double delta_h_ph = 0;
			double delta_p = 0;
			delta_p = delta_h*(last_ph+ph)/2.0;
			if( delta_h != 0) {
				double ph_h = delta_ph/delta_h;
				double ph0 = last_ph - ph_h*last_h;
				double int_0 = integral(ph0,ph_h,last_h);// ph_h*(last_h*last_h)/2.0 + ph0*last_h;
				double int_1 = integral(ph0,ph_h,h);//ph_h*(h*h)/2.0 + ph0*h;
				delta_h_ph = int_1 - int_0; 
				
				if( delta_h_ph == delta_h_ph) {
					total_h_ph += delta_h_ph;
				}
				
				if( ph_h != 0) {
					double delta_e = integral_e_of_e(ph0,ph_h,h) - integral_e_of_e(ph0,ph_h,last_h);
					if( delta_e != delta_e) {
						/*
						System.out.println("nan: "+ph0+","+ph_h+","+h+","+last_h);
						System.out.println(" : "+(ph0+ph_h*h));
						System.out.println(" : "+(ph0+ph_h*last_h));
						*/
						
					} else {
						total_e -= delta_e;
					}
				}
	
			}
			last_ph = ph;
			last_h = h;
		}
		//total_e /= //Math.log(2);
		return new double[]{total_h_ph,0-total_e,total_ph};
	}
	public static double integral(double a, double b, double x) {
		return x*x*(2.0*b*x+3.0*a)/6.0;
	}
	
	public static double integral_e_of_e(double a, double b, double x) {
		//integral of (a+bx)log(a+bx)
		double y = a+b*x;
		if( y == 0) {
			y = Double.MIN_VALUE;
		}
		return y*y*(2.0*Math.log(y)-1.0) / (4.0*b);
	}
	public Vector<double[]> bin(Vector<double[]> samples, int num_bins, int h_col, int p_col, double maxH) {
		int samples_per_bin = samples.size()/num_bins;
		Vector<double[]> bins = new Vector<double[]>();

		double lower = 0;
		for( int i = 0; i < num_bins; i++) {
			double p = 0;
			for( int j = 0; j < samples_per_bin; j++) {
				p += samples.get(i*samples_per_bin+j)[p_col];
			}

			double upper;
			if( i + 1 == num_bins ) {
				upper = maxH;
			} else {
				double uleft = samples.get((i+1)*samples_per_bin-1)[h_col];
				double uright = i + 1 == num_bins ? maxH : samples.get((i+1)*samples_per_bin)[h_col];
				upper = (uleft+uright)/2.0;
				
			}
			bins.add(new double[]{(lower+upper)/2.0,p});
			lower = upper;
		}
		return bins;
	}
	public Vector<double[]> getEntropyPDF(Integer[] data, int num_samples, int num_bins) {
		setNumberOfCategories(data.length);
		Vector<double[]> results = new Vector<double[]>();
		Vector<Pair<Double,double[]>> samples = new Vector<Pair<Double,double[]>>();
		
		//do a bunch of random samples
		for( int i = 0; i < num_samples; i++) {
			
			//calculate entropy and probability
			double[] dd = getRandomSample(data,false);
			samples.add(new Pair<Double,double[]>(dd[0],dd));
		}
		
		//now sort by entropy
		Collections.sort(samples);
		for(int i = 0; i < samples.size(); i++) {
			results.add(samples.get(i).b);
		}
		double maxH = Math.log((double)data.length);/////Math.log(2);

		Vector<double[]> dd = bin(results,num_bins,0,1,maxH);
		//Vector<double[]> samples, int num_bins, int h_col, int p_col) {
		
		return dd;
	}
	public double[] getRandomThetas(Integer[] ii) {
		if( !use_prior_for_sampling) {
			return getRandomThetas(ii.length, false);
		}
		double[] thetas = new double[ii.length];
		
		double sum = 0;
		for(int j = 0; j < thetas.length; j++) {
			thetas[j] = new org.apache.commons.math3.distribution.GammaDistribution((double)ii[j], 1.0).sample();
			sum += thetas[j];
		}
		
		//normalize so sum(p)=1
		if( sum == 0) {
			sum = 1;
		}
		sum = 1/sum;
		for(int j = 0; j < thetas.length; j++) {
			thetas[j] *= sum;
		}
		for(int j = 0; j < thetas.length; j++) {
			double r = thetas[j];
			int k = (int)(Math.random()*(double)thetas.length);
			thetas[j] = thetas[k];
			thetas[k] = r;
		}
		return thetas;
	}
	public double[] getRandomThetas(int n,boolean b) {
		double[] thetas = new double[n];
		double sum = 0;
		for(int j = 0; j < thetas.length; j++) {
			thetas[j] = 1.0/(Math.random());
			sum += thetas[j];
		}
		
		//normalize so sum(p)=1
		if( sum == 0) {
			sum = 1;
		}
		sum = 1/sum;
		for(int j = 0; j < thetas.length; j++) {
			thetas[j] *= sum;
		}
		return thetas;
	}
	public double[] getRandomSample(Integer[] data, boolean permuted) {
		return getEntropyAndProbabilityAtTheta(getRandomThetas(data),data,permuted);
	}
	public double[] getEntropyAndProbabilityAtTheta(double[] thetas, Integer[] data) {
		return getEntropyAndProbabilityAtTheta(thetas, data,true); 
	}
	public double[] getEntropyAndProbabilityAtTheta(double[] thetas, Integer[] data, boolean permuted) {
		
		double h = getEntropyGivenTheta(thetas);
		double pDO = 0;
		if( permuted) {
			pDO = getProbabilityOfDataGivenThetaPermuted2(thetas, data);
		} else {
			pDO = getProbabilityOfDataGivenTheta(thetas, data);
		}
		//double pDO = Math.exp(lpDO);
		
		
		double[] dHdOs = getDerivsOfEntropyGivenTheta(thetas);
		double dHdO = 1;//getDerivOfEntropyGivenTheta(thetas);
		for( int i = 0; i < dHdOs.length-1; i++) {
			dHdO *= dHdOs[i];
		}
		
		
		//double dHdO = determinant(getJacobianGivenTheta(thetas));
		//double lpO = Math.log(prior.getPriorProbabilityOfTheta(thetas));
		//double lph = lpDO+lpO - ( prior_is_on_H ? 0 : Math.log(Math.abs(dHdO)));
		
		//double l_add = sum_data;//data.length == 0 ? 0 : Math.log(data.length)*sum_data;
		//double l_add = data.length == 0 ? 0 : Math.log(data.length)*sum_data;
		//double l_add = data.length == 0 ? 0 : data.length*sum_data;
		//l_add -= 1000;
		//lph += l_add;
		
		//double pO = Math.exp(lpO);
		//double ph = Math.exp(lph);
		double pO = Math.log(prior.getPriorProbabilityOfTheta(thetas));
		double ph = pDO*pO / ( prior_is_on_H ? 1 : Math.abs(dHdO));
		//System.out.println("dHdO: "+dHdO+" pO: "+pO+" pDO: "+pDO+" ph: "+ph+" h: "+h+" p is on h: "+prior_is_on_H);
		/*
		if( prior_is_on_H)
			pO = prior.getPriorProbabilityOfTheta(h);
			*/
		if( ph == 0) {
			//System.out.println("ph is zero 25! "+pDO+" "+pO+" "+dHdO+" "+lpO+" "+lph+" "+l_add+" "+(lph+l_add));
		}
		return new double[]{h,ph,pDO,dHdO,pO};
	}
	public double[] getEntropyAndLogProbabilityAtTheta(double[] thetas, Integer[] data, boolean permuted) {

		
		double h = getEntropyGivenTheta(thetas);
		double lpDO = 0;
		if( permuted) {
			lpDO = Math.log(getProbabilityOfDataGivenThetaPermuted2(thetas, data));
		} else {
			lpDO = getLogProbabilityOfDataGivenTheta(thetas, data);
		}
		//double pDO = Math.exp(lpDO);
		
		
		double[] dHdOs = getDerivsOfEntropyGivenTheta(thetas);
		double dHdO = 1;//getDerivOfEntropyGivenTheta(thetas);
		for( int i = 0; i < dHdOs.length-1; i++) {
			dHdO *= dHdOs[i];
		}
		
		
		//double dHdO = determinant(getJacobianGivenTheta(thetas));
		double lpO = Math.log(prior.getPriorProbabilityOfTheta(thetas));
		double lph = lpDO+lpO - ( prior_is_on_H ? 0 : Math.log(Math.abs(dHdO)));
		
		//double l_add = sum_data;//data.length == 0 ? 0 : Math.log(data.length)*sum_data;
		//double l_add = data.length == 0 ? 0 : Math.log(data.length)*sum_data;
		//double l_add = data.length == 0 ? 0 : data.length*sum_data;
		//l_add -= 1000;
		//lph += l_add;
		
		//double pO = Math.exp(lpO);
		//double ph = Math.exp(lph);
		//System.out.println("dHdO: "+dHdO+" pO: "+pO+" pDO: "+pDO+" ph: "+ph+" h: "+h+" p is on h: "+prior_is_on_H);
		/*
		if( prior_is_on_H)
			pO = prior.getPriorProbabilityOfTheta(h);
			*/
		//if( ph == 0) {
			//System.out.println("ph is zero 25! "+pDO+" "+pO+" "+dHdO+" "+lpO+" "+lph+" "+l_add+" "+(lph+l_add));
		//}
		return new double[]{h,lph,lpDO,dHdO,lpO};
	}

	public double getDerivOfEntropyGivenTheta(double[] thetas) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double getEntropyGivenTheta(double[] theta) {
		double p = 0;
		double sumtheta = 0;
		for( int i = 0; i < theta.length; i++) {
			sumtheta += theta[i];
		}
		for( int i = 0; i < theta.length; i++) {
			double t = theta[i] / sumtheta;
			p -= t == 0 ? 0 : t*Math.log(t);
		}
		return p;

	}
	/*
	 * d(entropy)/d(1-a) = ln(1-a) - [b*ln(b) + c*ln(c)]/a
	 * 
	 * 
	 * d(entropy)/d(1-a) = ln(1-a) - [b*ln(b) + c*ln(c)]/a - a = 0 : ln(1)=0 - ...[current entropy = ?] /0

	 * d(entropy)/d(1-a) = ln(1-a) - [b*ln(b) + c*ln(c)]/a - a = 1 : ln(0)=inf - ...[current entropy = 0] /1
	 * 
	 * where a = b+c+....
	 * 
	 */
	
	// http://mathworld.wolfram.com/ChangeofVariablesTheorem.html
	
    public static double determinant(double A[][]) {
    	return determinant(A,A.length);
    }
    public static double determinant(double A[][],int N) {
        double det=0;
        if(N == 1) {
            det = A[0][0];
        } else 
        if (N == 2) {
            det = A[0][0]*A[1][1] - A[1][0]*A[0][1];
        } else {
            det=0;
            for(int j1=0;j1<N;j1++) {
                double[][] m = new double[N-1][];
                for(int k=0;k<(N-1);k++) {
                    m[k] = new double[N-1];
                }
                for(int i=1;i<N;i++) {
                    int j2=0;
                    for(int j=0;j<N;j++) {
                        if(j == j1) {
                            continue;
                        }
                        m[i-1][j2] = A[i][j];
                        j2++;
                    }
                }
                det += Math.pow(-1.0,1.0+j1+1.0)* A[0][j1] * determinant(m,N-1);
            }
        }
        return det;
    }
	 

	public double[][] getJacobianGivenTheta(double[] thetas) {
		boolean print = true;
		double[] partial_entropies = new double[thetas.length];
		double[][] jacob = new double[thetas.length][thetas.length];
		
		for( int i = 0; i < partial_entropies.length; i++) {
			partial_entropies[i] = (thetas[i] == 0 || thetas[i] == 1) ? 0 : -thetas[i] * Math.log(thetas[i]); 
		}
		for( int i = 0; i < thetas.length; i++) {
			if( thetas[i] == 1) {
				for( int j = 0; j < thetas.length; j++) {
					jacob[i][j] = (i == j) ? 0 : Double.POSITIVE_INFINITY;
				}
			} else
			if( thetas[i] == 0) {
				for( int j = 0; j < thetas.length; j++) {
					jacob[i][j] = (i == j) ? Double.NEGATIVE_INFINITY : 0;
				}
			} else {
				double a = 1-thetas[i];
				for( int j = 0; j < thetas.length; j++) {
					jacob[i][j] = (i == j) ? -Math.log(thetas[i]) : -partial_entropies[j]/a;
				}
			}
		}
		return jacob;
	}

	public double[] getDerivsOfEntropyGivenTheta(double[] thetas) {
		boolean print = false;
		double[] derivs = new double[thetas.length];
		if( print) System.out.print("[");
		
		for( int i = 0; i < derivs.length; i++) {
			if( print) System.out.print((i > 0 ? "," : "")+thetas[i]);
		}
		if( print) System.out.print("] : [");
		
		double base_entropy = getEntropyGivenTheta(thetas);
		
		for( int i = 0; i < derivs.length; i++) {
			if( thetas[i] == 0) {
				derivs[i] = Double.NEGATIVE_INFINITY;
				continue;
			} else
			if( thetas[i] == 1) {
				derivs[i] = Double.POSITIVE_INFINITY;
				continue;
			} 
			double partial_entropy = - thetas[i] * Math.log(thetas[i]);
			
			double right_part = (base_entropy - partial_entropy)/(1-thetas[i]);
			double left_part = Math.log(thetas[i]);
			
			derivs[i] = -(left_part+right_part);//thetas[i] == 0 ? 1 : -Math.log(thetas[i])-1;
			if( print) System.out.print((i > 0 ? "," : "")+derivs[i]);
		}
		if( print) System.out.println("]");
		return derivs;
	}
	/*//old way
	public double[] getDerivsOfEntropyGivenTheta(double[] thetas) {
		double[] derivs = new double[thetas.length];
		//System.out.print("[");
		for( int i = 0; i < derivs.length; i++) {
			//System.out.print((i > 0 ? "," : "")+thetas[i]);
		}
		//System.out.print("] : [");
		for( int i = 0; i < derivs.length; i++) {
			derivs[i] = thetas[i] == 0 ? 1 : -Math.log(thetas[i])-1;
			//System.out.print((i > 0 ? "," : "")+derivs[i]);
		}
		//System.out.println("]");
		return derivs;
	}
	*/

	@Override
	public double getProbabilityOfDataGivenTheta(double[] thetas, Integer[] data) {
		double mult = 1;
		for(int i = 0; i < thetas.length; i++) {
			mult *= Math.pow(thetas[i],(double)data[i]);
		}
		return mult;
		// TODO Auto-generated method stub
		//return Math.exp(getLogProbabilityOfDataGivenTheta(thetas,data));
	}

	public double getProbabilityOfThetaGivenData(double[] theta, Integer[] data) {
		return getProbabilityOfDataGivenTheta(theta,data)*prior.getPriorProbabilityOfTheta(theta);
	}
		
	public double getLogProbabilityOfDataGivenTheta(double[] thetas, Integer[] data) {
		/*
		double mult = 1;
		for(int i = 0; i < thetas.length; i++) {
			mult *= Math.pow(thetas[i],(double)data[i]);
		}
		return mult;
		*/
		
		
		double sum = 0;
		//double tot = 0;
		for(int i = 0; i < thetas.length; i++) {
			sum += Math.log(thetas[i])*(double)data[i];
			//tot += data[i];
		}
		//sum /= tot;
		return sum;
		//return Math.exp(sum);
		
		
	}

	public PriorDistribution<double[], double[], Integer> getConjugatePrior() {
		return new DirichletDistribution();
	}

	public double getProbabilityOfDataGivenThetaPermuted(double[] theta, Integer[] data) {
		if( permutations == null) {
			generatePermutations(theta.length);
		}
		double p = 0;
		double[] permuted = new double[theta.length];
		for( int i = 0; i < permutations.length; i++) {
			for( int j = 0; j < theta.length; j++) {
				permuted[j] = theta[permutations[i][j]];
			}
			p += getProbabilityOfDataGivenTheta(permuted,data); 
		}
		return p;
	}
	
	public double getProbabilityOfDataGivenThetaPermuted2(double[] theta, Integer[] data) {
		int bucket_count = theta.length;
		int[] buckets = new int[bucket_count];
		for( int i = 0; i < bucket_count; i++) {
			buckets[i] = i;
		}
		permutations = new int[(int)factorial(bucket_count)][];
		//Double p = new Double(0);
		permute2(buckets, 0, 0,data, new Double(0),theta);
		
		
		double p = 0;
		double[] permuted = new double[theta.length];
		for( int i = 0; i < permutations.length; i++) {
			for( int j = 0; j < theta.length; j++) {
				permuted[j] = theta[permutations[i][j]];
			}
			p += getProbabilityOfDataGivenTheta(permuted,data); 
		}
		return p;
	}
	int permute2(int[] num, int i, int start, Integer[] data, Double p, double[] theta) {
		if (start >= num.length) {
			double[] result = new double[num.length];
			for( int j = 0; j < num.length; j++) {
				result[j] = theta[num[j]];
			}
			p += getProbabilityOfDataGivenTheta(result,data); 
			i++;
		}
	 
		for (int j = start; j <= num.length - 1; j++) {
			swap(num, start, j);
			i = permute2(num, i, start + 1,data,p,theta);
			swap(num, start, j);
		}
		return i;
	}

	void generatePermutations(int bucket_count) {
		int[] buckets = new int[bucket_count];
		for( int i = 0; i < bucket_count; i++) {
			buckets[i] = i;
		}
		permutations = new int[(int)factorial(bucket_count)][];
		permute(buckets, permutations, 0, 0);
	}
	static int permute(int[] num, int[][] result, int i, int start) {
	 
		if (start >= num.length) {
			result[i] = new int[num.length];
			for( int j = 0; j < num.length; j++)
				result[i][j] = num[j];
			i++;
		}
	 
		for (int j = start; j <= num.length - 1; j++) {
			swap(num, start, j);
			i = permute(num, result, i, start + 1);
			swap(num, start, j);
		}
		return i;
	}
	static void swap(int[] a, int i, int j) {
		int temp = a[i];
		a[i] = a[j];
		a[j] = temp;
	}

	static long factorial(long n) {
		if( n <= 1) {
			return 1;
		}
		return n*factorial(n-1);
	}
	
	public static double[] getMaxEntropyThetas(int n) {
		if( n == 0) {
			return new double[]{};
		}
		double p = 1.0/(double)n;
		double[] dd = new double[n];
		for( int i = 0; i < dd.length; i++) {
			dd[i] = p;
		}
		return dd;
	}
	public static double[] getMinEntropyThetas(int n) {
		if( n == 0) {
			return new double[]{};
		}
		//double p = 1.0/(thetas.length);
		double[] dd = new double[n];
		for( int i = 0; i < dd.length; i++) {
			dd[i] = 0;
		}
		dd[0] = 1;
		return dd;
	}
	
	public int getRandomBayesian(Integer[] buckets) {
		return 0;
	}
	public int getRandomFrequentist(Integer[] buckets) {
		return 0;
	}
	
	//also known as marginal likelihood
	//https://en.wikipedia.org/wiki/Categorical_distribution
	//�Dirichlet-Multinomial� 
	//https://people.eecs.berkeley.edu/~stephentu/writeups/dirichlet-conjugate-prior.pdf
	public double[] getPDFBayesian(Integer[] buckets, int num_samples) {
		Integer[][] data = new Integer[buckets.length][];
		for( int i = 0; i < data.length; i++) {
			data[i] = new Integer[data.length];
			data[i][i] = 1;
		}
		double[] dd = new double[thetas.length];
		for( int i = 0; i < num_samples; i++) {
			double[] t = getRandomThetas(thetas.length,false);
			for( int j = 0; j < dd.length; j++) {
				dd[j] += getProbabilityOfThetaGivenData(t,buckets) * getProbabilityOfDataGivenTheta(t,data[j]);
			}
		}
		
		double sum = 0;
		for( int j = 0; j < dd.length; j++) {
			sum += dd[j];
		}
		sum = 1.0/sum;
		for( int j = 0; j < dd.length; j++) {
			dd[j] *= sum;
		}
		return dd;
	}




}
