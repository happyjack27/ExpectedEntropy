package util;

import org.apache.commons.math3.special.*;
import org.apache.commons.math3.util.*;

public class Functions {
	
	//  D_Paisley_John_a_201005.pdf
	
	/*************************************
	 An ANSI-C implementation of the digamma-function for real arguments based
	 on the Chebyshev expansion proposed in appendix E of 
	 http://arXiv.org/abs/FastMath.CA/0403344 . This is identical to the implementation
	 by Jet Wimp, FastMath. Comp. vol 15 no 74 (1961) pp 174 (see Table 1).
	 For other implementations see
	 the GSL implementation for Psi(Digamma) in
	 http://www.gnu.org/software/gsl/manual/html_node/Psi-_0028Digamma_0029-Function.html

	Richard J. Mathar, 2005-11-24
	**************************************/
	//#include <FastMath.h>

	//#ifndef M_PIl
	/** The constant Pi in high precision */
	//#define M_PIl 3.1415926535897932384626433832795029L
	//#endif

	/** Euler's constant in high precision */
	public static final double M_GAMMAl = 0.5772156649015328606065120900824024;
	/** the natural logarithm of 2 in high precision */
	public static final double  M_LN2l = 0.6931471805599453094172321214581766;

	/** The digamma function in long double precision.
	* @param x the real value of the argument
	* @return the value of the digamma (psi) function at that point
	* @author Richard J. Mathar
	* @since 2005-11-24
	*/

	
	
	public static double GAMMA = 0.577215664901532860606512090082;
	public static double GAMMA_MINX = 1.e-12;
	public static double DIGAMMA_MINNEGX = -1250;
	public static double C_LIMIT = 49;
	public static double S_LIMIT = 1e-5;
	
	public static double digammal( double x)
	{
		/* force into the interval 1..3 */
		if( x < 0.0 )
			return digammal(1.0-x)+FastMath.PI/FastMath.tan(FastMath.PI*(1.0-x)) ;	/* reflection formula */
		else if( x < 1.0 )
			return digammal(1.0+x)-1.0/x ;
		else if ( x == 1.0)
			return -M_GAMMAl ;
		else if ( x == 2.0)
			return 1.0-M_GAMMAl ;
		else if ( x == 3.0)
			return 1.5-M_GAMMAl ;
		else if ( x > 3.0)
			/* duplication formula */
			return 0.5*(digammal(x/2.0)+digammal((x+1.0)/2.0))+M_LN2l ;
		else
		{
			/* Just for your information, the following lines contain
			* the Maple source code to re-generate the table that is
			* eventually becoming the Kncoe[] array below
			* interface(prettyprint=0) :
			* Digits := 63 :
			* r := 0 :
			* 
			* for l from 1 to 60 do
			* 	d := binomial(-1/2,l) :
			* 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
			* 	evalf(r) ;
			* 	print(%,evalf(1+Psi(1)-r)) ;
			*o d :
			* 
			* for N from 1 to 28 do
			* 	r := 0 :
			* 	n := N-1 :
			*
	 		*	for l from iquo(n+3,2) to 70 do
			*		d := 0 :
	 		*		for s from 0 to n+1 do
	 		*		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
	 		*		od :
	 		*		if 2*l-n > 1 then
	 		*		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
	 		*		fi :
	 		*	od :
	 		*	print(evalf((-1)^n*2*r)) ;
	 		*od :
	 		*quit :
			*/
			double Kncoe[] = { .30459198558715155634315638246624251,
			.72037977439182833573548891941219706, -.12454959243861367729528855995001087,
			.27769457331927827002810119567456810e-1, -.67762371439822456447373550186163070e-2,
			.17238755142247705209823876688592170e-2, -.44817699064252933515310345718960928e-3,
			.11793660000155572716272710617753373e-3, -.31253894280980134452125172274246963e-4,
			.83173997012173283398932708991137488e-5, -.22191427643780045431149221890172210e-5,
			.59302266729329346291029599913617915e-6, -.15863051191470655433559920279603632e-6,
			.42459203983193603241777510648681429e-7, -.11369129616951114238848106591780146e-7,
			.304502217295931698401459168423403510e-8, -.81568455080753152802915013641723686e-9,
			.21852324749975455125936715817306383e-9, -.58546491441689515680751900276454407e-10,
			.15686348450871204869813586459513648e-10, -.42029496273143231373796179302482033e-11,
			.11261435719264907097227520956710754e-11, -.30174353636860279765375177200637590e-12,
			.80850955256389526647406571868193768e-13, -.21663779809421233144009565199997351e-13,
			.58047634271339391495076374966835526e-14, -.15553767189204733561108869588173845e-14,
			.41676108598040807753707828039353330e-15, -.11167065064221317094734023242188463e-15 } ;

			double Tn_1 = 1.0 ;	/* T_{n-1}(x), started at n=1 */
			double Tn = x-2.0 ;	/* T_{n}(x) , started at n=1 */
			double resul = Kncoe[0] + Kncoe[1]*Tn ;

			x -= 2.0 ;

			for(int n = 2 ; n < Kncoe.length ;n++)
			{
				double Tn1 = 2.0 * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
				resul += Kncoe[n]*Tn1 ;
				Tn_1 = Tn ;
				Tn = Tn1 ;
			}
			return resul ;
		}
	}


public static double digammalAlt( double x, int n)
{
	/* force into the interval 1..3 */
	if( x < 0.0 )
		return digammalAlt(1.0-x,n)+FastMath.PI/FastMath.tan(FastMath.PI*(1.0-x)) ;	/* reflection formula */
	else if( x < 1.0 )
		return digammalAlt(1.0+x,n)-1.0/x ;
	else if ( x == 1.0)
		return -M_GAMMAl ;
	else if ( x == 2.0)
		return 1.0-M_GAMMAl ;
	else if ( x == 3.0)
		return 1.5-M_GAMMAl ;
	else if ( x > 3.0)
		return digammalAlt(x-1.0,n)+1.0/(x-1.0) ;
	else
	{
		x -= 1.0;
		double resul = -M_GAMMAl ;

		for( ; n >= 1 ;n--)
			resul += x/(n*(n+x)) ;
		return resul ;
	}
}

/*
public static double getExpectedEntropyForDelta(double cur_expected_entropy, double a, double g_count, double g_count_delta) {  //remember to add +1 before passing in
	double ra = 1.0/a;
	double a_new = a + g_count_delta;
	double ra_new = 1/a_new;
	double g = g_count*ra;
	double g_new = (g_count+g_count_delta)*ra_new;
	
	double old_result_part = digamma(a+1)-g*digamma(a*g+1);
	double new_result_part = digamma(a_new+1)-g_new*digamma(a_new*g_new+1);

	return cur_expected_entropy-old_result_part+new_result_part;
}*/


//this uses a prior on p, rather than entropy!  i.e. entropy(1,1)=0.72... instead of 0.5
public static double getExpectedEntropy(Integer[] integers, double exaggerate, double prior) {
	double a = 0;
	double[] p = new double[integers.length];
	for( int i = 0; i < p.length; i++) {
		a += (integers[i]*exaggerate+prior);
	}
	double ra = 1.0/a;
	for( int i = 0; i < p.length; i++) {
		p[i] = ra*(double)(integers[i]*exaggerate+prior);
	}
	double result = digamma(a+1);
	
	for( int i = 0; i < p.length; i++) {
		result -= p[i]*digamma((integers[i]*exaggerate+prior)+1);
	}
	return result;
}
	
public static double getNaiveExpectedEntropy(double num_buckets, double prior) {
	return digamma(num_buckets*prior+1) - digamma(prior+1);
}
	
public static double getExpectedEntropyVariance(Integer[] buckets) {
	double a = 0;
	double[] g = new double[buckets.length];
	double[] ag = new double[buckets.length];
	for( int i = 0; i < g.length; i++) {
		a += (buckets[i]+1);
	}
	double ra = 1.0/a;
	for( int i = 0; i < g.length; i++) {
		g[i] = ra*(double)(buckets[i]+1);
	}
	for( int i = 0; i < g.length; i++) {
		ag[i] = buckets[i]+1;
	}
	
	double sum1 = 0;
	for( int i = 0; i < g.length; i++) {
		double amt = (
				FastMath.pow(digamma(ag[i]+2)-digamma(a+2), 2)
				+Gamma.trigamma(ag[i]+2)-Gamma.trigamma(a+2)
				)
				*g[i]*(ag[i]+1)/(a+1)
				;
		sum1 += amt;//g[i]*Functions.digamma((buckets[i]+1)+1);
	}
	
	double sum2 = 0;
	for( int i = 0; i < g.length; i++) {
		for( int j = 0; j < g.length; j++) {
			if( i == j) {
				continue;
			}
			double amt = (
					(digamma(ag[i]+1)-digamma(a+2))
					*(digamma(ag[j]+1)-digamma(a+2)-Gamma.trigamma(a+2))
					)
					*a*g[i]*g[j]/(a+1)
					;
			sum2 += amt;//g[i]*Functions.digamma((buckets[i]+1)+1);
		}
	}
	
	double expected_entropy = getExpectedEntropy(buckets,1,1);
	return sum1+sum2-expected_entropy*expected_entropy;
}
	

	public static double logGamma(double x) {
		double tmp = (x - 0.5) * FastMath.log(x + 4.5) - (x + 4.5);
		double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
		                       + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
		                       +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
		return tmp + FastMath.log(ser * FastMath.sqrt(2 * FastMath.PI));
	}
	public static double gamma(double x) { return FastMath.exp(logGamma(x)); }
	public static double digamma(double x) {

	    double value = 0;

	    while (true){

	        if (x >= 0 && x < GAMMA_MINX) {
	            x = GAMMA_MINX;
	        }
	        if (x < DIGAMMA_MINNEGX) {
	            x = DIGAMMA_MINNEGX + GAMMA_MINX;
	            continue;
	        }
	        if (x > 0 && x <= S_LIMIT) {
	            return value + -GAMMA - 1 / x;
	        }

	        if (x >= C_LIMIT) {
	            double inv = 1 / (x * x);
	            return value + FastMath.log(x) - 0.5 / x - inv
	                    * ((1.0 / 12) + inv * (1.0 / 120 - inv / 252));
	        }

	        value -= 1 / x;
	        x = x + 1;
	    }

	}	
}
