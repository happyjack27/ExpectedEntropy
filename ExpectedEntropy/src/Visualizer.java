
public class Visualizer {
	/*
	 * 
	 * stage 1: make a graph of connections
	 * 
	 * calculate H(w,x,y,z), that's how many e's per dot.
	 * 
	 * calculate Dwxyz = I(w,x,y,z), connect w,x,y,z to that. those dots are used up.
	 * calculate Dwxy  = I(w,x,y)-Dwxyz, connect w,x,y to that. those dots are used up.
	 * ...
	 * calculate Dwx = I(w,x)-Dwxy-Dwxz-Dwxyz, connect w,x to that. those dots are used up.
	 * ...
	 * calculate Dw = I(w) - all higher cardinality sets that include w, connect w to that. those dots are used up.
	 * 
	 */
	/*
	 * 
	 * stage 2: sort the graph
	 * 
	 * the sum of the squared distance between each category's dots squared distance from their centroid is the "fitness score" of the map.
	 * 
	 * use genetic crossover, with mutation performed by randomly perturbing every point's position according to a normal distribution, with standard deviation decaying exponentially over time.
	 * 
	 * 	 * or, with a grid system, flip a neighboring pair with probability p.  or choose n pairs at random to flip where n is choosen from a poisson distribution, with exponentially decaying rate parameter. 
	 * 
	 * 
	 */
	
	int num_vars = 4;
	int num_dots = 100*100;
	int last_dot = 0;
	double totalH = 4;
	double dots_per_H = 0;
	double[] I;
	double[] D;
	int size = 0;
	int[] dot_connections = new int[num_dots];
	int[][] dot_coords = new int[num_dots][2];
	
	public void init() {
		size = 0x01 << num_vars;
		I = new double[size];
		D = new double[size];
		dots_per_H = (double)num_dots/totalH;
		
		//calculate dot_density
		for( int bit_count = num_vars; bit_count > 0; bit_count--) {
			for( int i = 0; i < size; i++) {
				if( Integer.bitCount(i) == bit_count) {
					D[i] = I[i] - totalParents(i);
				}
			}
		}

		//now attach dots
		int i = 0;
		double dot_remainder = 0;
		double cur_remainder = D[i]*dots_per_H;
		for( int d = 0; d < num_dots; d++) {
			if( cur_remainder < 1) {
				dot_remainder += cur_remainder;
				if( dot_remainder >= 1) {
					dot_connections[d] = i;
					dot_remainder--;
				} else {
					d--;
				}
				i++;
				cur_remainder = D[i]*dots_per_H;
			} else {
				dot_connections[d] = i;
				cur_remainder--;
			}
		}
	}
	public double scoreMap(int[][] dot_coords) {
		double[] counts = new double[num_vars];
		double[][] centers = new double[num_vars][2];
		double[] squared_distances = new double[num_vars];
		
		//compute centers
		for( int d = 0; d < dot_connections.length; d++) {
			for( int i = 0; i < num_vars; i++) {
				if( ((0x01 << i) & dot_connections[d]) != 0 ) {
					counts[i]++;
					centers[i][0] += dot_coords[i][0];
					centers[i][1] += dot_coords[i][1];
				}
			}
		}
		for( int i = 0; i < num_vars; i++) {
			counts[i] = 1.0 / counts[i];
			centers[i][0] *= counts[i];
			centers[i][1] *= counts[i];
		}
		
		//compute sum squared distance
		double ssd = 0;
		for( int d = 0; d < dot_connections.length; d++) {
			for( int i = 0; i < num_vars; i++) {
				if( ((0x01 << i) & dot_connections[d]) != 0 ) {
					double dx = dot_coords[i][0] - centers[i][0];
					double dy = dot_coords[i][1] - centers[i][1];
					ssd += dx*dx + dy*dy;
				}
			}
		}
		return ssd;
	}
	
	public double totalParents(int mask) {
		double tot = 0;
		for( int i = 0; i < size; i++) {
			if( (i & mask) == mask && i != mask) {
				tot += D[i];
			}
		}
		return tot;
	}

}
