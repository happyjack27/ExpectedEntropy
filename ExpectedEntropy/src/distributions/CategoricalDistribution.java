package distributions;

import java.util.Vector;

import distributions.interfaces.PosteriorDistribution;
import distributions.interfaces.PriorDistribution;
import util.*;

public class CategoricalDistribution implements PosteriorDistribution<double[],double[],Integer> {
	double[] thetas;
	public boolean prior_is_on_H = false;
	int[][] permutations;
	
	PriorDistribution<double[],double[],Integer> prior;
	
	public static double GAMMA = 0.577215664901532860606512090082;
	public static double GAMMA_MINX = 1.e-12;
	public static double DIGAMMA_MINNEGX = -1250;
	public static double C_LIMIT = 49;
	public static double S_LIMIT = 1e-5;
	
	public static int[] bucket(Vector<boolean[]> samples, int[] select) {
		int[] buckets = new int[(int)Math.pow(2, select.length)];
		for( boolean[] s : samples) {
			int bucket = 0;
			for( int i = 0; i < select.length; i++) {
				bucket = bucket * 2 + (s[select[i]] == true ? 1 : 0);
			}
			buckets[bucket]++;
		}
		return buckets;
	}
	static int[][] sets = new int[][]{
		new int[]{0},
		new int[]{1},
		new int[]{2},
		new int[]{3},
		new int[]{1,2,3},//5
		new int[]{0,2,3},
		new int[]{0,1,3},//7
		new int[]{0,1,2},
		new int[]{1,2,3,0},
		new int[]{0,3},
		new int[]{1,3},
		new int[]{2,3},
	};
	
	public static void main(String[] args) {
		int[] buckets = new int[16];
		Vector<boolean[]> samples = new Vector<boolean[]>();
		
		for( int i = 0; i < 256; i++) {
			boolean a = Math.random() > 0.5;
			boolean b = Math.random() > 0.5;
			boolean c = a ^ b;
			boolean d = a & b;
			samples.add(new boolean[]{a,b,c,d});
			//System.out.println(a+" "+b+" "+c+" "+d+" ");
			
			for(int[] ss : sets) {
				double v = Functions.getExpectedEntropy(bucket(samples,ss),1,1)/Math.log(2);
				System.out.print( v+", ");
			}
			System.out.println();
		}
		
		/*
		System.out.println("{0,0}: "+Functions.getExpectedEntropy(new int[]{0,0},1,1)/Math.log(2));
		System.out.println();
		System.out.println("{0,1}: "+Functions.getExpectedEntropy(new int[]{0,1},1,1)/Math.log(2));
		System.out.println();
		System.out.println("{0,2}: "+Functions.getExpectedEntropy(new int[]{0,2},1,1)/Math.log(2));
		System.out.println("{1,1}: "+Functions.getExpectedEntropy(new int[]{1,1},1,1)/Math.log(2));
		System.out.println();
		System.out.println("{0,3}: "+Functions.getExpectedEntropy(new int[]{0,3},1,1)/Math.log(2));
		System.out.println("{1,2}: "+Functions.getExpectedEntropy(new int[]{1,2},1,1)/Math.log(2));
		System.out.println();
		System.out.println("{0,4}: "+Functions.getExpectedEntropy(new int[]{0,4},1,1)/Math.log(2));
		System.out.println("{1,3}: "+Functions.getExpectedEntropy(new int[]{1,3},1,1)/Math.log(2));
		System.out.println("{2,2}: "+Functions.getExpectedEntropy(new int[]{2,2},1,1)/Math.log(2));

		System.out.println();
		System.out.println("{200,200}: "+Functions.getExpectedEntropy(new int[]{200,200},1,1)/Math.log(2));
		System.out.println("{0,400}: "+Functions.getExpectedEntropy(new int[]{0,400},1,1)/Math.log(2));
		*/
		
	}
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

	public double[] getEntropyAndProbabilityAtTheta(double[] thetas, Integer[] data) {
		double h = getEntropyGivenTheta(thetas);
		double pDO = getProbabilityOfDataGivenThetaPermuted(thetas, data);
		
		
		double[] dHdOs = getDerivsOfEntropyGivenTheta(thetas);
		double dHdO = 1;//getDerivOfEntropyGivenTheta(thetas);
		for( int i = 0; i < dHdOs.length-1; i++) {
			dHdO *= dHdOs[i];
		}
		
		
		//double dHdO = determinant(getJacobianGivenTheta(thetas));
		double pO = prior.getPriorProbabilityOfTheta(thetas);
		double ph = pDO*pO / ( prior_is_on_H ? 1 : Math.abs(dHdO));
		System.out.println("dHdO: "+dHdO+" pO: "+pO+" pDO: "+pDO+" ph: "+ph+" h: "+h+" p is on h: "+prior_is_on_H);
		/*
		if( prior_is_on_H)
			pO = prior.getPriorProbabilityOfTheta(h);
			*/
		
		return new double[]{h,ph,pDO,dHdO,pO};
	}

	public double getDerivOfEntropyGivenTheta(double[] thetas) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double getEntropyGivenTheta(double[] theta) {
		double p = 0;
		for( int i = 0; i < theta.length; i++) {
			p -= theta[i] == 0 ? 0 : theta[i]*Math.log(theta[i]);
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
		boolean print = true;
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

	public double getProbabilityOfDataGivenTheta(double[] thetas, Integer[] data) {
		double mult = 1;
		for(int i = 0; i < thetas.length; i++) {
			mult *= Math.pow(thetas[i],(double)data[i]);
		}
		return mult;
		/*
		double sum = 0;
		for(int i = 0; i < thetas.length; i++) {
			sum += Math.log(thetas[i])*(double)data[i];
		}
		return Math.exp(sum);
		*/
	}

	public double getProbabilityOfThetaGivenData(double[] theta, Integer[] data) {
		return getProbabilityOfDataGivenTheta(theta,data)*prior.getPriorProbabilityOfTheta(theta);
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
	


}
