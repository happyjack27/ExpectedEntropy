
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
	int dot_width = 64;
	int num_dots = dot_width*dot_width;
	int last_dot = 0;
	double totalH = 4;
	double dots_per_H = 0;
	double[] I;
	double[] D;
	int size = 0;
	int[] dot_connections = new int[num_dots];
	int[][] dot_coords = new int[num_dots][2];
	
	public void init(int num, double totalEntropy, double[] Is) {
		num_vars = num;
		totalH = totalEntropy;
		size = 0x01 << num_vars;
		I = Is;
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
		
		//now randomize array
		for( int d = 0; d < dot_connections.length; d++) {
			int d2 = (int)(Math.random()*(double)dot_connections.length);
			int t = dot_connections[d];
			dot_connections[d] = dot_connections[d2];
			dot_connections[d2] = t;
		}
		
	}
	
	public void perturb(int[][] dot_grid, double rate) {
		int N = (int)rate; //better to poisson estimate this.
		for( int i = 0; i < N; i++) {
			int x = (int)(Math.random()*(double)dot_grid.length);
			int y = (int)(Math.random()*(double)dot_grid.length);
			int dx = (int)(Math.random()*3.0-1.0); 
			int dy = (int)(Math.random()*3.0-1.0); 
			int x2 = (x+dx) < 0 ? x : (x+dx) >= dot_grid.length ? x : (x+dx);
			int y2 = (y+dy) < 0 ? y : (y+dy) >= dot_grid.length ? y : (y+dy);
			int t = dot_grid[x][y];
			dot_grid[x][y] = dot_grid[x2][y2];
			dot_grid[x2][y2] = t;
		}
	}
	
	public double scoreGrid(int[][] dot_grid) {
		double[] counts = new double[num_vars];
		double[][] centers = new double[num_vars][2];
		double[] squared_distances = new double[num_vars];
		
		//compute centers
		for( int x = 0; x < dot_grid.length; x++) {
			for( int y = 0; y < dot_grid.length; y++) {
				int dot = dot_grid[x][y];
				for( int i = 0; i < num_vars; i++) {
					if( ((0x01 << i) & dot) != 0 ) {
						counts[i]++;
						centers[i][0] += x;//dot_coords[i][0];
						centers[i][1] += y;//dot_coords[i][1];
					}
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
		for( int x = 0; x < dot_grid.length; x++) {
			for( int y = 0; y < dot_grid.length; y++) {
				int dot = dot_grid[x][y];
				for( int i = 0; i < num_vars; i++) {
					if( ((0x01 << i) & dot) != 0 ) {
						double dx = x - centers[i][0];
						double dy = y - centers[i][1];
						ssd += dx*dx + dy*dy;
					}
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
