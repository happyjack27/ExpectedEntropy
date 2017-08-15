
import java.util.*;

import distributions.CategoricalDistribution;
import util.*;

// http://www.sortie-nd.org/lme/Statistical%20Papers/Burnham_and_Anderson_2004_Multimodel_Inference.pdf

public class _RunTests {
	static int max_n = 25;
	static int num_runs = 50;//1000;//25;
	static boolean divide_log_likelihood_by_n = false;
	static int MONTE_CARLO_RESOLUTION = 100;
	
	static int number_of_bits = 10;
	static double number_of_nats = number_of_bits*Math.log(2.0);
	public static int max_size = 10;

	static double at_percentile = 0;
	static boolean show_detail = false;
	
	
	//entropy at percentile
	//percentile at entropy
	//expected entropy
	
	
	
	static final int METRIC_BEES = 0;
	static final int METRIC_BEES_PRIOR_NOT_H = 1;
	static final int METRIC_BEES_START_WITH_1 = 8;
	static final int METRIC_AIC = 2;
	static final int METRIC_LOG_LIKELIHOOD = 9;
	static final int METRIC_BIC = 4;

	static final int METRIC_AICc = 3;
	static final int METRIC_BEES_LNMULT = 5;
	static final int METRIC_BEES_I_MULT = 6;
	static final int METRIC_BEES_C_MULT = 7;
	static final int METRIC_BEES_PENALIZED_LOG = 11;
	static final int METRIC_BEES_PENALIZED_K_N = 12;
	static final int METRIC_BEES_PENALIZED_K_N_NEG = 13;
	static final int METRIC_BEES_PENALIZED_MULT = 14;
	static final int METRIC_BEES_CONSTRAINED = 15;
	//static double min_entropy_reduction_per_parameter = 0.1;
	static double min_entropy_reduction_per_parameter_for_k_n = 1.0;
	//static double min_e_per_p = 0.1;
	static double constrained_min_total_e_saved_per_p = 1;
	
	static boolean use_neg_entropy = true;
	
	static boolean use_prior_on_p = false;

	
	static double C = 10;
	 
	static int METRIC = 0;


	//0 - wire
	//1 - and/or/xor
	//2 - half add
	//8 - full add
	public static int[] choices = new int[]{
			0,1,2,8,
			//0,1,2,8,
			//0,1,2,8,
			//0,1,2,8,
			
			//0,1,2,8,
			//0,1,2,8,
			//0,1,2,8,
			//8,8,8,
			//8,8,8,
			//2,2,
			//2,2,2,
			//2,2,2,2,
			//0,1,0,1,
			//3,3,3,3,3,
			/*
			0,0,0,0,0,
			1,1,1,1,1,
			0,0,0,0,0,
			1,1,1,1,1,
			*/
			//0,0,0,0,
			//1,1,1,1,
			/*
			0,0,0,0,//0,0,0,0,
			1,1,1,1,//1,1,1,1,

			0,0,0,0,//0,0,0,0,
			1,1,1,1,//1,1,1,1,

			0,0,0,0,//0,0,0,0,
			1,1,1,1,//1,1,1,1,  

			0,0,0,0,//0,0,0,0,
			1,1,1,1,//1,1,1,1,
			*/
	};
	public static double[] noises = new double[]{
			/*
			0.1,0.1,0.1,0.1,
			0.1,0.1,0.1,0.1,
			0.1,0.1,0.1,0.1,
			0.1,0.1,0.1,0.1,
			0.1,0.1,0.1,0.1,
			*/
			
			0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,
			0.1,0.1,0.1,0.1,
			//0.0,0.0,0.1,0.1,
			/*
			0.0,0.0,0.0,
			0.05,0.05,0.1,0.1,0.2,0.2,
			*/
			//0.1,0.1,
			
			
			//0.0,0.0,0.0, 0.0,0.0,0.0,
			0.1,0.1,0.1, 0.1,0.1,0.1,
			0.0,0.0,0.0,0.0,0.0,0.0,
			//0.1,0.1,0.1,0.1,0.1,0.1,
			
			
			//0.05,0.05,0.05,0.05,//0.0,0.0,0.0,0.0,
			//0.05,0.05,0.05,0.05,//0.0,0.0,0.0,0.0,

			0.1,0.1,0.1,0.1,//0.1,0.1,0.1,0.1,
			0.1,0.1,0.1,0.1,//0.1,0.1,0.1,0.1,

			0.2,0.2,0.2,0.2,//0.2,0.2,0.2,0.2,
			0.2,0.2,0.2,0.2,//0.2,0.2,0.2,0.2,
	};
	public static int[] metrics = new int[]{
			//METRIC_BEES_START_WITH_1,METRIC_BEES_START_WITH_1,METRIC_BEES_START_WITH_1,METRIC_BEES_START_WITH_1,
			//METRIC_AIC,
			//METRIC_BEES,METRIC_BEES,METRIC_BEES,METRIC_BEES,
			//METRIC_BEES,METRIC_BEES,METRIC_BEES,METRIC_BEES,
			//METRIC_LOG_LIKELIHOOD,METRIC_LOG_LIKELIHOOD,METRIC_LOG_LIKELIHOOD,METRIC_LOG_LIKELIHOOD,
			//METRIC_AICc,METRIC_AICc,METRIC_AICc,METRIC_AICc,
			//METRIC_BIC,METRIC_BIC,METRIC_BIC,METRIC_BIC,
			
			METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,

			METRIC_LOG_LIKELIHOOD,METRIC_LOG_LIKELIHOOD,METRIC_LOG_LIKELIHOOD,METRIC_LOG_LIKELIHOOD,
			METRIC_BEES,METRIC_BEES,METRIC_BEES,METRIC_BEES,
			METRIC_AIC,METRIC_AIC,METRIC_AIC,METRIC_AIC,
			METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,
			//METRIC_BEES_PENALIZED_MULT,METRIC_BEES_PENALIZED_MULT,METRIC_BEES_PENALIZED_MULT,METRIC_BEES_PENALIZED_MULT,
			//METRIC_BEES,
			//METRIC_LOG_LIKELIHOOD,
			//METRIC_BEES_PRIOR_NOT_H,METRIC_BEES_PRIOR_NOT_H,METRIC_BEES_PRIOR_NOT_H,METRIC_BEES_PRIOR_NOT_H,
			//METRIC_LOG_LIKELIHOOD,METRIC_AIC,
			//METRIC_BEES,METRIC_BEES,METRIC_BEES,
			METRIC_LOG_LIKELIHOOD,METRIC_AIC,METRIC_BEES,
			METRIC_LOG_LIKELIHOOD,METRIC_AIC,METRIC_BEES,
			METRIC_BEES,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_MULT,METRIC_BEES_PRIOR_NOT_H,
			//METRIC_LOG_LIKELIHOOD,METRIC_AIC,METRIC_BEES, METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_MULT,METRIC_BEES_PRIOR_NOT_H,

			//METRIC_LOG_LIKELIHOOD,METRIC_AIC, 
			METRIC_BEES,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_MULT,METRIC_BEES_PRIOR_NOT_H,
			//METRIC_LOG_LIKELIHOOD,METRIC_AIC,METRIC_BEES, METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_MULT,METRIC_BEES_PRIOR_NOT_H,

			/*
			METRIC_BEES_PENALIZED_MULT,METRIC_BEES_CONSTRAINED,
			METRIC_LOG_LIKELIHOOD,METRIC_AIC,METRIC_BEES,METRIC_BEES_PENALIZED_K_N,
			METRIC_LOG_LIKELIHOOD,METRIC_AIC,METRIC_BEES,METRIC_BEES_PENALIZED_K_N,
			*/

			//METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,
			//METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PENALIZED_K_N,
			//METRIC_BEES_PENALIZED_LOG,METRIC_BEES_PENALIZED_LOG,
			//METRIC_BEES_START_WITH_1,METRIC_BEES_PENALIZED_LOG,
			/*
			METRIC_BEES,METRIC_BEES_PENALIZED,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PRIOR_NOT_H,
			METRIC_BEES,METRIC_BEES_PENALIZED,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PRIOR_NOT_H,
			METRIC_BEES,METRIC_BEES_PENALIZED,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PRIOR_NOT_H,
			METRIC_BEES,METRIC_BEES_PENALIZED,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PRIOR_NOT_H,
			METRIC_BEES,METRIC_BEES_PENALIZED,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PRIOR_NOT_H,
			METRIC_BEES,METRIC_BEES_PENALIZED,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PRIOR_NOT_H,
			METRIC_BEES,METRIC_BEES_PENALIZED,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PRIOR_NOT_H,
			METRIC_BEES,METRIC_BEES_PENALIZED,METRIC_BEES_PENALIZED_K_N,METRIC_BEES_PRIOR_NOT_H,
			*/
			/*
			METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_AICc,METRIC_BIC,
			METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_AICc,METRIC_BIC,
			METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_AICc,METRIC_BIC,
			METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_AICc,METRIC_BIC,
			METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_AICc,METRIC_BIC,
			METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_AICc,METRIC_BIC,
			METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_AICc,METRIC_BIC,
			METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_AICc,METRIC_BIC,
			*/
			//METRIC_BEES,METRIC_BEES_PENALIZED_LOG,METRIC_BEES_PENALIZED,METRIC_BEES_START_WITH_1,METRIC_BEES_PRIOR_NOT_H,METRIC_AIC_NO_PENALTY,METRIC_AIC,METRIC_BIC,
			
	};
	
	public static void setTargetCover() {
		target_cover = new Vector<Integer>();
		if( max_size == 3) {
			switch(choice) {
			case 0:
				target_cover.add(3);
				target_cover.add(12);
				target_cover.add(37);
				target_cover.add(86);
				target_cover.add(167);
				break;
			case 1:
				
				target_cover.add(7);
				target_cover.add(41);
				target_cover.add(129);
				target_cover.add(130);
				break;
			case 2:
				/*
				target_cover.add(6);
				target_cover.add(56);
				target_cover.add(448);
				target_cover.add(513);
				*/
				//sample[9];
				break;
			case 3:
				/*
				target_cover.add(1);
				target_cover.add(30);
				target_cover.add(992);
				*/
				break;
			}		
		}
		if( max_size == 4) {
			switch(choice) {
			case 0:
				break;
			case 1:
				break;
			case 2:
				target_cover.add(30);
				target_cover.add(385);
				target_cover.add(1);
				target_cover.add(31);
				//sample[9];
				break;
			case 3:
				break;
			}		
		}
		if( max_size == 5) {
			switch(choice) {
			case 0:
				break;
			case 1:
				break;
			case 2:
				//sample[9];
				break;
			case 3:
				break;
			case 8:
				target_cover.add(31);
				target_cover.add(637);
				break;
			}		
		}
		if( max_size < 10) {
			return;
		}
		switch(choice) {
		case 0:
			target_cover.add(3);
			target_cover.add(12);
			target_cover.add(48);
			target_cover.add(192);
			target_cover.add(768);
			break;
		case 1:
			target_cover.add(7);
			target_cover.add(56);
			target_cover.add(448);
			target_cover.add(512);
			break;
		case 2:
			
			target_cover.add(30);
			target_cover.add(960);
			target_cover.add(1);
			target_cover.add(32);
			//sample[9];
			break;
			/*
		case 3:
			target_cover.add(1);
			target_cover.add(30);
			target_cover.add(992);
			break;
			*/
		case 8:
			target_cover.add(31);
			target_cover.add(992);
			//target_cover.add(1);
			//target_cover.add(30);
			//target_cover.add(992);
			break;
		}		
	}
	
	
	
	public static int cur_choice = 0;
	public static int choice = 0;
	public static double noise = 0.0;

	public static boolean optimize_cover = true;
	static Vector<Integer> target_cover = new Vector<Integer>();

	static int[][] sets = null;	
	static Vector<Vector<Integer>> all_covers = new Vector<Vector<Integer>>();	
	
	static Vector<double[]> runs = new Vector<double[]>(); 
	static double[] results = new double[max_n];
	static double[] run_avg = new double[max_n];
	static double[][] run_avgs = new double[choices.length][max_n];
	
	public static Vector<Pair<Double,Vector<Integer>>> covers = new Vector<Pair<Double,Vector<Integer>>>();
	
	public static void main( String[] args) {
		for( cur_choice = 0; cur_choice < choices.length; cur_choice++) {
			/*
			at_percentile = 
					cur_choice == 0 ? 0.25 :
					cur_choice == 1 ? 0.50 :
					cur_choice == 2 ? 0.75 :
					0
			;
			*/
			
			choice = choices[cur_choice];
			noise = noises[cur_choice];
			METRIC = metrics[cur_choice];
			
			Vector<int[]> vsets = allSubSets(number_of_bits,max_size);
			sets = new int[vsets.size()][];
			for(int i = 0; i < sets.length; i++) {
				sets[i] = vsets.get(i);
			}
			addAllCovers();
			System.out.println("found "+all_covers.size()+" unique covers");
			
			
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
				
			for( int i = 0; i < num_runs; i++) {
				doRun();
				runs.add(results);
				results = new double[max_n];
			}
			double m = 1.0/(double)num_runs;
			for( int i = 0; i < max_n; i++) {
				run_avg[i]*=m;
				System.out.println(run_avg[i]);
			}
			run_avgs[cur_choice] = run_avg.clone();
			System.out.println();
		}
		System.out.println();
		for( int i = 0; i < max_n; i++) {
			for( int j = 0; j < run_avgs.length; j++) {
				System.out.print(run_avgs[j][i]+", ");
			}
			System.out.println();
		}
	}
	
	public static void addAllCovers() {
		all_covers = new Vector<Vector<Integer>>();
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
	
	public static Vector<Integer> scoreAndStuff( double[] entropies, double N) {
		double best_e = 999999999999999999999999999.9;
		double best_e_param_count = 1;
		Vector<Integer> best_cover = new Vector<Integer>();
		/*
	static final int METRIC_BEES = 0;
	static final int METRIC_BEES_PRIOR_NOT_H = 1;
	static final int METRIC_BEES_START_WITH_1 = 8;
		 */
		for( Vector<Integer> cover : all_covers) {
			double e = getTotalEntropy(cover,entropies);
			double p = getTotalParams(cover);
			double max_e = best_e;
			//if( use_neg_entropy || METRIC == METRIC_BEES_PENALIZED_K_N_NEG) {
				e -= number_of_nats;
			//}
			e *= N;
			
			
			if( METRIC == METRIC_BEES_CONSTRAINED) {
				if( -e / p < constrained_min_total_e_saved_per_p) {
					e /= p;
				}
			}
			if( METRIC == METRIC_BEES_PENALIZED_MULT) {
				e /= p;//Math.log(p)*min_entropy_reduction_per_parameter_for_k_n;
			}
			if( METRIC == METRIC_BEES_PENALIZED_LOG) {
				e += Math.log(p)*min_entropy_reduction_per_parameter_for_k_n;
			}
			if( METRIC == METRIC_BEES_PENALIZED_K_N) {
				//if( METRIC == METRIC_BEES || METRIC == METRIC_BEES_PRIOR_NOT_H || METRIC == METRIC_BEES_START_WITH_1) {
				//max_e -= (p-best_e_param_count)*min_entropy_reduction_per_parameter_for_k_n / N;
				e += p*min_entropy_reduction_per_parameter_for_k_n;
				//e += min_entropy_reduction_per_parameter_for_k_n*p / N;
			}
			//min_entropy_reduction_per_parameter
			if( e < max_e || e == max_e && p < best_e_param_count) {
				best_e = e;
				best_e_param_count = p;
				best_cover = cover;
			}
		}
		return best_cover;

	}
	
	private static double getTotalParams(Vector<Integer> cover) {
		double tot = 0;
		for( int i : cover) {
			tot += 0x01 << sets[i].length;
		}
		return tot;
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
		boolean[] sample = new boolean[number_of_bits];
		for( int i = 0; i < sample.length; i++) {
			sample[i] = Math.random() > 0.5;
		}
		
		switch(choice) {
		case 0:
			sample[1] = sample[0];
			sample[3] = sample[2];
			sample[5] = sample[4];
			sample[7] = sample[6];
			sample[9] = sample[8];
			break;
		case 1:
			sample[2] = sample[0] ^ sample[1];
			sample[5] = sample[4] & sample[3];
			sample[8] = sample[6] | sample[7];
			break;
		case 2:
			//sample[0]; //mixed //2 half adders
			sample[3] = sample[1] ^ sample[2];
			sample[4] = sample[1] & sample[2];

			sample[8] = sample[6] ^ sample[7];
			sample[9] = sample[6] & sample[7];
			//sample[9];
			break;
		case 3:
			//half adder, full adder
			sample[3] = sample[1] ^ sample[2];
			sample[4] = sample[1] & sample[2];

			sample[8] = sample[5] ^ sample[6] ^ sample[7];
			sample[9] = (sample[5] & sample[6]) | (sample[6] & sample[7]) | (sample[7] & sample[5]);
			break;
		case 4:
			//something random
			sample[2] = sample[0] ^ sample[1];
			sample[3] = sample[0] & sample[1];
			sample[6] = sample[5];
			sample[9] = sample[7] | sample[8];
			break;
		case 8:
			sample[3] = sample[0] ^ sample[1] ^ sample[2];
			sample[4] = (sample[0] & sample[1]) | (sample[1] & sample[2]) | (sample[2] & sample[0]);
			
			sample[8] = sample[5] ^ sample[6] ^ sample[7];
			sample[9] = (sample[5] & sample[6]) | (sample[6] & sample[7]) | (sample[7] & sample[5]);
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
	
	public static double getAICNoPenalty(Integer[] buckets) {
		double k = buckets.length;
		
		double total_log_probability = 0;
		
		double n = 0;
		for(int i = 0; i < buckets.length; i++) {
			n += buckets[i];
		}
		
		for(int i = 0; i < buckets.length; i++) {
			double p = ((double)buckets[i]) / n;
			if( p == 0) {
				continue;
			}
			total_log_probability += Math.log(p)*(double)buckets[i];
		}
		if( divide_log_likelihood_by_n) {
			total_log_probability /= n;
		}
		
		return  - total_log_probability;
	}
	
	public static double getAIC(Integer[] buckets) {
		double k = buckets.length;
		
		double total_log_probability = 0;
		
		double n = 0;
		for(int i = 0; i < buckets.length; i++) {
			n += buckets[i];
		}
		
		for(int i = 0; i < buckets.length; i++) {
			double p = ((double)buckets[i]) / n;
			if( p == 0) {
				continue;
			}
			total_log_probability += Math.log(p)*(double)buckets[i];
		}
		if( divide_log_likelihood_by_n) {
			total_log_probability /= n;
		}
		
		return k - total_log_probability;
	}
	public static double getAICc(Integer[] buckets) {
		double k = buckets.length;
		
		double total_log_probability = 0;
		
		double n = 0;
		for(int i = 0; i < buckets.length; i++) {
			n += buckets[i];
		}
		for(int i = 0; i < buckets.length; i++) {
			double p = ((double)buckets[i]) / n;
			if( p == 0) {
				continue;
			}
			total_log_probability += Math.log(p)*(double)buckets[i];
		}
		if( divide_log_likelihood_by_n) {
			total_log_probability /= n;
		}
		
		if( n-k-1.0 == 0) {
			return k - total_log_probability;
		}
		return k - total_log_probability + k*(k+1.0)/(n-k-1.0);
	}
	
	public static double getBIC(Integer[] buckets) {
		double k = buckets.length;
		
		double total_log_probability = 0;
		
		double n = 0;
		for(int i = 0; i < buckets.length; i++) {
			n += buckets[i];
		}
		for(int i = 0; i < buckets.length; i++) {
			double p = ((double)buckets[i]) / n;
			if( p == 0) {
				continue;
			}
			total_log_probability += Math.log(p)*(double)buckets[i];
		}
		
		if( divide_log_likelihood_by_n) {
			total_log_probability /= n;
		}
		
		return Math.log((double)n)*k - 2.0*total_log_probability;
	}
	

	public static void doRun() {
		//System.exit(0);
		setTargetCover();
		Vector<Integer> best_cover = getRandomFullCover();
		CategoricalDistribution cat = new CategoricalDistribution();
		
		int num_bins = 10240;
		int samples_per_bin = 10240;
		int[] buckets = new int[16];
		Vector<boolean[]> samples = new Vector<boolean[]>();
		
		for( int i = 0; i < max_n; i++) {
			samples.add(getSample());
			
			double[] entropies = new double[sets.length];
			for(int j = 0; j < sets.length; j++) {//
				int[] ss = sets[j];
				Integer[] ii = bucket(samples,ss);
				cat.setNumberOfCategories(ii.length);
				double v = 0;
				cat.prior_is_on_H = true;
				
				switch(METRIC) {
				case METRIC_BEES:
				case METRIC_BEES_PENALIZED_K_N:
				case METRIC_BEES_PENALIZED_LOG:
				case METRIC_BEES_PENALIZED_MULT:
				case METRIC_BEES_CONSTRAINED:
				case METRIC_BEES_PENALIZED_K_N_NEG:
					if( use_prior_on_p) {
						v = Functions.getExpectedEntropy(ii,1,1);///Math.log(2);
					} else {
						Vector<double[]> curve = cat.getEntropyCurveMultiplicative(ii,MONTE_CARLO_RESOLUTION);
						if( at_percentile > 0) {
							v = cat.getAtPercentile(curve, at_percentile);
						} else {
							v = cat.integrateSortedEntropyCurve(curve);
						}
					}

					
					//CategoricalDistribution cat = new CategoricalDistribution();
					//v = cat.getSummaryStats(ii,MONTE_CARLO_RESOLUTION)[0];
					//System.out.print(""+ii.length+": "+v);
					//System.out.print(".");
					//v = Functions.getExpectedEntropy(ii,1,1)/Math.log(2);
					break;
				case METRIC_BEES_C_MULT:
					Integer[] ii1 = new Integer[ii.length];
					double mult1 = C;
					for( int k = 0; k < ii.length; k++) {
						ii1[k] = (int)(mult1*(double)ii[k]);
					}
					v = cat.getSummaryStats(ii1,MONTE_CARLO_RESOLUTION)[0];
					break;
				case METRIC_BEES_I_MULT:
					Integer[] ii2 = new Integer[ii.length];
					double mult = ii.length;
					for( int k = 0; k < ii.length; k++) {
						ii2[k] = (int)(mult*(double)ii[k]);
					}
					v = cat.getSummaryStats(ii2,MONTE_CARLO_RESOLUTION)[0];
					break;
				case METRIC_BEES_LNMULT:
					Integer[] ii3 = new Integer[ii.length];
					double mult3 = Math.log((double)ii.length)/Math.log(2.0);
					for( int k = 0; k < ii.length; k++) {
						ii3[k] = (int)(mult3*(double)ii[k]);
					}
					v = cat.getSummaryStats(ii3,MONTE_CARLO_RESOLUTION)[0];
					//System.out.print(".");
					//v = Functions.getExpectedEntropy(ii,1,1)/Math.log(2);
					break;
				case METRIC_BEES_START_WITH_1:
					Integer[] ii4 = new Integer[ii.length];
					//double mult3 = Math.log((double)ii.length)/Math.log(2.0);
					for( int k = 0; k < ii.length; k++) {
						ii4[k] = ii[k]+1;
					}
					if( use_prior_on_p) {
						v = Functions.getExpectedEntropy(ii4,1,1);///Math.log(2);
					} else {
						Vector<double[]> curve = cat.getEntropyCurveMultiplicative(ii4,MONTE_CARLO_RESOLUTION);
						if( at_percentile > 0) {
							v = cat.getAtPercentile(curve, at_percentile);
						} else {
							v = cat.integrateSortedEntropyCurve(curve);
						}
					}
					break;
				case METRIC_BEES_PRIOR_NOT_H:
					v = Functions.getExpectedEntropy(ii,1,1);///Math.log(2);
					break;
				case METRIC_AIC:
					v = getAIC(ii);
					break;
				case METRIC_AICc:
					v = getAICc(ii);
					break;
				case METRIC_BIC:
					v = getBIC(ii);
					break;
				case METRIC_LOG_LIKELIHOOD:
					v = getAICNoPenalty(ii);
					break;
				}
				
				entropies[j] = v;
			}
			
			best_cover = scoreAndStuff(entropies,samples.size());
			double best_e = getTotalEntropy(best_cover,entropies);

			

			
			int num = getDistance(best_cover,target_cover);
			if( show_detail) {
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
				System.out.print(": "+num);
				System.out.println();
			}
			
			

			if( i % 100 == 0) {
				System.out.print(".");
			}
			results[i] = num;
			run_avg[i] += num;

			if( num == 0) {
				//System.exit(0);
			}
		}

		
	}
}
