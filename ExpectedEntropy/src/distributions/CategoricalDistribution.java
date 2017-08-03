package distributions;

import java.util.*;

import distributions.interfaces.PosteriorDistribution;
import distributions.interfaces.PriorDistribution;
import util.*;

public class CategoricalDistribution implements PosteriorDistribution<double[],double[],Integer> {

	public static int choice = 3;
	public static double noise = 0.10;

	double[] thetas;
	public boolean prior_is_on_H = true;
	public static boolean optimize_cover = true;
	static Vector<Integer> target_cover = new Vector<Integer>();
	static int number_of_bits = 10;
	public static int max_size = 10;

	int[][] permutations;
	static int[][] sets = null;	
	static Vector<Vector<Integer>> all_covers = new Vector<Vector<Integer>>();	
	
	PriorDistribution<double[],double[],Integer> prior;
	
	public static double GAMMA = 0.577215664901532860606512090082;
	public static double GAMMA_MINX = 1.e-12;
	public static double DIGAMMA_MINNEGX = -1250;
	public static double C_LIMIT = 49;
	public static double S_LIMIT = 1e-5;
	
	public static Vector<Pair<Double,Vector<Integer>>> covers = new Vector<Pair<Double,Vector<Integer>>>();
	
	public static void addAllCovers() {
		addAllCovers(new Vector<Integer>(), 0, 1); //skip the empty set
	}
	public static void addAllCovers(Vector<Integer> current, int covered, int next_set) {
		int full_cover = (0x01 << number_of_bits)-1;
		for( int set = next_set; set < sets.length; set++) {
			int[] ii = sets[set];
			int set_cover = 0;
			for( int i : ii) {
				set_cover |= 0x01 << i;
			}
			if( (set_cover & covered) != 0) {
				continue;
			}
			current.add(set);
			int new_cover = set_cover | covered;
			if( new_cover == full_cover) {
				Vector<Integer> full = (Vector<Integer>)current.clone();
				Collections.sort(full);;
				all_covers.add(full);
			} else {
				addAllCovers(current,new_cover,set+1);
			}
			current.remove(current.size()-1);
		}
	}
	
	public static Vector<Integer> scoreAndStuff( double[] entropies) {
		double best_e = 999999999999999999999999999.9;
		Vector<Integer> best_cover = new Vector<Integer>();
		for( Vector<Integer> cover : all_covers) {
			double e = getTotalEntropy(cover,entropies);
			if( e < best_e) {
				best_e = e;
				best_cover = cover;
			}
		}
		return best_cover;
		/*
		
		
		int population = 10000;
		
		//punish empty set
		entropies[0] = 1000;
		for( int n = 0; n < 5; n++) {
			while( covers.size() < population) {
				covers.add(new Pair<Double,Vector<Integer>>(0.0,getRandomFullCover()));
			}
			for( int i = 0; i < population; i++) {
				Vector<Integer> cover = (Vector<Integer>) covers.get(i).b.clone();
				cover.remove((int)(Math.random()*(double)cover.size()));
				while( !doesCover(cover)) {
					cover.add((int)(Math.random()*(double)sets.length));
				}
				if( optimize_cover) {
					optimizeFullCover(cover);
				}
				covers.add(new Pair<Double,Vector<Integer>>(0.0,cover));
			}
			for( Pair<Double,Vector<Integer>> c : covers) {
				c.a = getTotalEntropy(c.b,entropies);
			}
			Collections.sort(covers);
			while( covers.size() > population) {
				covers.remove(population);
			}
		}
		*/
	}
	
	public static Integer[] bucket(Vector<boolean[]> samples, int[] select) {
		if( select.length == 0) {
			return new Integer[0];
		}
		int pow = (int)(Math.pow(2, select.length));
		Integer[] buckets = new Integer[pow];
		for(int i = 0; i < buckets.length; i++) {
			buckets[i] = new Integer(0);
		}
		for( boolean[] s : samples) {
			int bucket = 0;
			for( int i = 0; i < select.length; i++) {
				bucket = bucket * 2 + (s[select[i]] == true ? 1 : 0);
			}
			buckets[bucket]++;
		}
		return buckets;
	}
	public static Vector<int[]> allSubSets(int n, int max_size) {
		Vector<int[]> sets = new Vector<int[]>();
		int combos = 0x01 << n;
		
		for( int i = 0; i < combos; i++) {
			if( Integer.bitCount(i) > max_size) {
				continue;
			}
			int[] a = new int[Integer.bitCount(i)];
			int index = 0;
			for( int j = 0; j < n; j++) {
				if( (i & (0x01 << j)) > 0) {
					a[index++] = j;
					//System.out.print(j+" ");
				}
			}
			//System.out.println();
			sets.add(a);
		}
		return sets;
	};
	public static int getDistance(Vector<Integer> cover, Vector<Integer> target_cover) {
		int d = 0;
		
		for(int s : cover) {
			int[] ii = sets[s];
			for( int i : ii) {
				int[] jj = sets[findSetWithMember(i,target_cover)];
				d += countMismatch(ii,jj,number_of_bits);
			}
			
		}
		
		return d;
	}
	public static int findSetWithMember(int m, Vector<Integer> cover) {
		for(int s : cover) {
			int[] ii = sets[s];
			for( int i : ii) {
				if( i == m) {
					return s;
				}
			}
			
		}
		return 0;
	}
	public static int countMismatch(int[] ii, int[] jj, int n) {
		int mm = 0;
		for( int i = 0; i < n; i++) {
			boolean ii_contains = false;
			for( int j = 0; j < ii.length; j++) {
				if( ii[j] == i) {
					ii_contains = true;
					break;
				}
			}
			boolean jj_contains = false;
			for( int j = 0; j < jj.length; j++) {
				if( jj[j] == i) {
					jj_contains = true;
					break;
				}
			}
			if( ii_contains ^ jj_contains) {
				mm++;
			}
		}
		return mm;
	}
	
	public static Vector<Integer> getRandomFullCover() {
		Vector<Integer> cover = new Vector<Integer>();
		while( !doesCover(cover)) {
			cover.add((int)(Math.random()*(double)sets.length));
		}
		optimizeFullCover(cover);
		return cover;
	}
	public static void optimizeFullCover(Vector<Integer> v) {
		//randomize first
		for( int i = 0; i < v.size(); i++) {
			int r = v.remove((int)(Math.random()*(double)v.size()));
			v.add(r);
		}
		for( int i = 0; i < v.size(); i++) {
			int r = v.remove(i);
			if( doesCover(v) ) {
				i--;
			} else {
				v.add(i, r);
			}
		}
		//now sort
		Collections.sort(v);
	}
	public static boolean doesCover(Vector<Integer> model) {
		boolean[] covered = new boolean[number_of_bits];
		for( int i = 0; i < covered.length; i++) {
			covered[i] = false;
		}
		for( Integer s : model) {
			for( int i : sets[s]) {
				covered[i] = true;
			}
		}
		
		boolean all_covered = true;
		for( int i = 0; i < covered.length; i++) {
			all_covered &= covered[i];
		}
		return all_covered;
	}
	/*
	static int[][] sets = new int[][]{
		new int[]{0},
		new int[]{1},
		new int[]{2},
		new int[]{3}, //sb 0.80
		new int[]{1,2,3},//sb 2
		new int[]{0,2,3},
		new int[]{0,1,3},
		new int[]{0,1,2},
		new int[]{1,2,3,0}, //sb 2
		new int[]{0,3}, //sb 1.5
		new int[]{1,3},
		new int[]{2,3},
	};
	*/
	
	
	public static boolean[] getSample() {
		boolean[] sample = new boolean[]{
				Math.random() > 0.5,
				Math.random() > 0.5,
				Math.random() > 0.5,
				Math.random() > 0.5,
				Math.random() > 0.5,
				Math.random() > 0.5,
				Math.random() > 0.5,
				Math.random() > 0.5,
				Math.random() > 0.5,
				Math.random() > 0.5,
		};
		
		switch(choice) {
		case 0:
			sample[1] = sample[0];
			sample[3] = sample[2];
			sample[5] = sample[4];
			sample[7] = sample[6];
			sample[9] = sample[8];
			target_cover.add(3);
			target_cover.add(12);
			target_cover.add(48);
			target_cover.add(192);
			target_cover.add(768);
			break;
		case 1:
			sample[2] = sample[0] ^ sample[1];
			sample[5] = sample[4] & sample[3];
			sample[8] = sample[6] | sample[7];
			target_cover.add(7);
			target_cover.add(56);
			target_cover.add(448);
			target_cover.add(512);
			break;
		case 2:
			//sample[0];
			sample[2] = sample[1];
			sample[5] = sample[4] & sample[3];
			sample[8] = sample[6] ^ sample[7];
			sample[9] = sample[0];
			target_cover.add(6);
			target_cover.add(56);
			target_cover.add(448);
			target_cover.add(513);
			//sample[9];
			break;
		case 3:
			//half adder, full adder
			sample[3] = sample[1] ^ sample[2];
			sample[4] = sample[1] & sample[2];

			sample[8] = sample[5] ^ sample[6] ^ sample[7];
			sample[9] = (sample[5] & sample[6]) | (sample[6] & sample[7]) | (sample[7] & sample[5]);
			target_cover.add(1);
			target_cover.add(30);
			target_cover.add(992);
			break;
		case 4:
			//something random
			sample[2] = sample[0] ^ sample[1];
			sample[3] = sample[0] & sample[1];
			sample[6] = sample[5];
			sample[9] = sample[7] | sample[8];
			break;
		case 8:
			//random
			break;
		case 9:
			//all 0s
			for(int i = 0; i < number_of_bits; i++) {
				sample[i] = false;
			}
			break;
		case 10:
			//all joined
			for(int i = 0; i < number_of_bits; i++) {
				sample[i] = sample[0];
			}
			break;
		}
		for( int i = 0; i < sample.length; i++) {
			if( Math.random() < noise) {
				sample[i] = Math.random() > 0.5;
			}
		}
		

		return sample;
		/*
		boolean a = Math.random() > 0.5;
		boolean b = Math.random() > 0.5;
		boolean c = a ^ b;
		boolean d = a & b;
		return new boolean[]{a,b,c,d};
		*/
	}
	public static double getTotalEntropy(Vector<Integer> cover, double[] entropies) {
		double tot = 0;
		for( int i : cover) {
			tot += entropies[i];
		}
		return tot;
	}
	
	public static void main(String[] args) {
		Vector<int[]> vsets = allSubSets(number_of_bits,max_size);
		sets = new int[vsets.size()][];
		for(int i = 0; i < sets.length; i++) {
			sets[i] = vsets.get(i);
		}
		addAllCovers();
		System.out.println("found "+all_covers.size()+" unique covers");
		
		
		//System.exit(0);
		Vector<Integer> best_cover = getRandomFullCover();
		
		int num_bins = 10240;
		int samples_per_bin = 10240;
		int[] buckets = new int[16];
		Vector<boolean[]> samples = new Vector<boolean[]>();
		
		CategoricalDistribution cat = new CategoricalDistribution();
		//cat.prior_is_on_H = true;
		/*
		System.out.println(""+cat.getEntropyGivenTheta(new double[]{1,0})/Math.log(2));
		System.out.println(""+cat.getEntropyGivenTheta(new double[]{1,1,0,0})/Math.log(2));
		System.out.println(""+cat.getEntropyGivenTheta(new double[]{1,1,0,0,0,0,0})/Math.log(2));
		System.out.println(""+cat.getEntropyGivenTheta(new double[]{1,1,1,1,0,0,0,0,0,0})/Math.log(2));
		*/
		//System.exit(0);
		
		//cat.setNumberOfCategories(8);
		/*
		int mult = 0;
		int o = 0;
		Integer[] ii2 = new Integer[]{o,o,mult,mult
				,o,o,o,o
				//,o,o,o,o,o,o,o,o
				};
		cat.setNumberOfCategories(ii2.length);
		Vector<double[]> vdd = cat.getEntropyPDF(ii2,num_bins*samples_per_bin, num_bins);
		double ddp = 0;
		for( int i = 0; i < vdd.size(); i++) {
			double[] dd = vdd.get(i);
			ddp += dd[1];
			System.out.println(""+dd[0]/Math.log(2)+", "+dd[1]+", "+ddp);
		}
		System.exit(0);
		*/
		
		for( int i = 0; i < 500; i++) {
			samples.add(getSample());
			/*
			if( i != 0) {
				for( int j = 0; j < 10000; j++) {
					samples.add(getSample());
				}
			}
			if( i == 0) {
				for( int j = 0; j < 1000; j++) {
					//samples.add(getSample());
				}
			}
			*/
			//System.out.println(a+" "+b+" "+c+" "+d+" ");
			
			double[] entropies = new double[sets.length];
			for(int j = 0; j < sets.length; j++) {//
				int[] ss = sets[j];
				Integer[] ii = bucket(samples,ss);
				cat.setNumberOfCategories(ii.length);
				//double[] dd = getSummaryStatsForEntropyPDF(cat.getEntropyPDF(ii,num_bins*samples_per_bin, num_bins));
				//System.out.print( dd[0]+", ");
				double v = Functions.getExpectedEntropy(ii,1,1)/Math.log(2);
				entropies[j] = v;
				//System.out.print( v+", ");
			}
			
			best_cover = scoreAndStuff(entropies);
			double best_e = getTotalEntropy(best_cover,entropies);
			//double best_e = covers.get(0).a;
			//best_cover = covers.get(0).b;
					/*
			getTotalEntropy(best_cover,entropies);
			for( int j = 0; j < 1000000; j++) {
				Vector<Integer> test_cover = getRandomFullCover();
				double test_e = getTotalEntropy(test_cover,entropies);
				if( test_e < best_e) {
					best_e = test_e;
					best_cover = test_cover;
				}
			}
			*/
			

			
			System.out.print(i+": best cover e: "+best_e+" : ");
			for( int s : best_cover) {
				int[] ii = sets[s]; 
				System.out.print(s+":{");
				//System.out.print("{");
				for( int j = 0; j < ii.length; j++) {
					if( j > 0)
						System.out.print(",");
						System.out.print(ii[j]);
				}
				System.out.print("} ");
			}
			

			int num = isCorrect(best_cover);
			System.out.print(": "+num);
			System.out.println();

			if( num == 0) {
				System.exit(0);
			}
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
	public static int isCorrect(Vector<Integer> cover) {
		return getDistance(cover,target_cover);
		/*
		int c = 0;
		c += cover.contains(378) ? 1 : 0;
		c += cover.contains(16) ? 1 : 0;
		c += cover.contains(15) ? 1 : 0;
		c += cover.contains(83) ? 1 : 0;
		return c;
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
	
	public static double[] getSummaryStatsForEntropyPDF( Vector<double[]> samples) {

		double last_h = 0;
		double last_ph = 0;
		double total_h_ph = 0;
		double total_ph = 0;
		double total_e = 0;
		
		//compute normalizing constant
		for( int i = 0; i < samples.size(); i++) {
			double[] result = samples.get(i);
			double h = result[0]/Math.log(2);
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
			double h = result[0]/Math.log(2);
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
		total_e /= Math.log(2);
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
		Vector<double[]> results = new Vector<double[]>();
		Vector<Pair<Double,double[]>> samples = new Vector<Pair<Double,double[]>>();
		
		//do a bunch of random samples
		for( int i = 0; i < num_samples; i++) {
			/*
			double sum = 0;
			double[] thetas = new double[data.length];
			for(int j = 0; j < thetas.length; j++) {
				thetas[j] = Math.random();
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
			*/
			
			//calculate entropy and probability
			double[] dd = getRandomSample(data,false);
			samples.add(new Pair<Double,double[]>(dd[0],dd));
		}
		
		//now sort by entropy
		Collections.sort(samples);
		for(int i = 0; i < samples.size(); i++) {
			results.add(samples.get(i).b);
		}
		double maxH = Math.log((double)data.length);///Math.log(2.0);

		Vector<double[]> dd = bin(results,num_bins,0,1,maxH);
		//Vector<double[]> samples, int num_bins, int h_col, int p_col) {
		
		return dd;
	}
	public double[] getRandomThetas(int n) {
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
		return getEntropyAndProbabilityAtTheta(getRandomThetas(data.length),data,permuted);
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
		
		
		double[] dHdOs = getDerivsOfEntropyGivenTheta(thetas);
		double dHdO = 1;//getDerivOfEntropyGivenTheta(thetas);
		for( int i = 0; i < dHdOs.length-1; i++) {
			dHdO *= dHdOs[i];
		}
		
		
		//double dHdO = determinant(getJacobianGivenTheta(thetas));
		double pO = prior.getPriorProbabilityOfTheta(thetas);
		double ph = pDO*pO / ( prior_is_on_H ? 1 : Math.abs(dHdO));
		//System.out.println("dHdO: "+dHdO+" pO: "+pO+" pDO: "+pDO+" ph: "+ph+" h: "+h+" p is on h: "+prior_is_on_H);
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
	
	public double getProbabilityOfDataGivenThetaPermuted2(double[] theta, Integer[] data) {
		int bucket_count = theta.length;
		int[] buckets = new int[bucket_count];
		for( int i = 0; i < bucket_count; i++) {
			buckets[i] = i;
		}
		//permutations = new int[(int)factorial(bucket_count)][];
		Double p = new Double(0);
		permute2(buckets, 0, 0,data,p,theta);
		
		/*
		double p = 0;
		double[] permuted = new double[theta.length];
		for( int i = 0; i < permutations.length; i++) {
			for( int j = 0; j < theta.length; j++) {
				permuted[j] = theta[permutations[i][j]];
			}
			p += getProbabilityOfDataGivenTheta(permuted,data); 
		}*/
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
	
	public int getRandomBayesian(Integer[] buckets) {
		return 0;
	}
	public int getRandomFrequentist(Integer[] buckets) {
		return 0;
	}
	
	//also known as marginal likelihood
	//https://en.wikipedia.org/wiki/Categorical_distribution
	//“Dirichlet-Multinomial” 
	//https://people.eecs.berkeley.edu/~stephentu/writeups/dirichlet-conjugate-prior.pdf
	public double[] getPDFBayesian(Integer[] buckets, int num_samples) {
		Integer[][] data = new Integer[buckets.length][];
		for( int i = 0; i < data.length; i++) {
			data[i] = new Integer[data.length];
			data[i][i] = 1;
		}
		double[] dd = new double[thetas.length];
		for( int i = 0; i < num_samples; i++) {
			double[] t = getRandomThetas(thetas.length);
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
