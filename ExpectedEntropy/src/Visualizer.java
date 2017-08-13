import java.awt.*;
import java.io.*;
import java.util.*;

import util.Pair;

public class Visualizer implements Draws {
	
	double totalH = 4;
	double dots_per_H = 0;
	double[] I;
	double[] D;
	int size = 0;

	int num_vars = 4;
	int dot_width = 64;
	int num_dots = dot_width*dot_width;
	int last_dot = 0;
	int[] dot_connections = new int[num_dots];
	int[][] dot_coords = new int[num_dots][2];
	static int[][] dot_grid = null;
	static Vector<Pair<Double,int[][]>> all_grids = new Vector<Pair<Double,int[][]>>();
	int num_grids = 100;
	
	double init_rate = 2;
	double anneal_mult = 0.99;
	static int image_size = 512;
	
	public static final int OUTPUT_FILE = 0;
	public static final int INPUT_FILE = 1;
	public static final int VARS = 2;
	public static final int GRIDSIZE = 3;
	public static final int ITERATIONS = 4;
	public static final int IMAGE_SIZE = 5;
	
	public static boolean use_squared_distance = true;
	public static boolean divide_by_area = false;
	public static double whitespace_fraction = 0.75;
	public static boolean shrink_by_area = false;
	
	public static void main(String[] args) {
		if( args.length > 4) {
			try {
				int gridSize = Integer.parseInt(args[GRIDSIZE]);
				int iterations = Integer.parseInt(args[ITERATIONS]);
				int numVars = Integer.parseInt(args[VARS]);
				File f = new File(args[INPUT_FILE]);
				if( args.length > 5) {
					image_size = Integer.parseInt(args[IMAGE_SIZE]);
				}
				
				FileInputStream fis = new FileInputStream(f);
				BufferedReader buf = new BufferedReader(new InputStreamReader(fis));
				double[] Is = new double[0x01 << numVars];
				for( int i = 0; i < Is.length; i++) {
					Is[i] = Double.parseDouble(buf.readLine());
				}
				Visualizer v = new Visualizer();
				v.init(numVars,gridSize,iterations,Is);
				
				dot_grid = v.all_grids.get(0).b;
				v.toCsv(dot_grid, args[OUTPUT_FILE]);
				
				ToImageFile img = new ToImageFile();
				System.out.println("drawing...");
				img.toPNG("visualizer.png", new Visualizer(), image_size);
				System.out.println("done.");

			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.exit(0);
		}
		/*
		double[] Is = new double[]{
				1.5,
				1,
				1,
				0.5
		};*/
		
		double[] Is = new double[]{
				3.00,
				1,
				1,
				0.5,
				2,
				0.25,
				0.25,
				0,
				
		};
		
		Visualizer v = new Visualizer();
		v.init(3,640,640,Is);
		
		
		ToImageFile img = new ToImageFile();
		System.out.println("drawing...");
		img.toPNG("visualizer.png", v, 640);
		System.out.println("done.");
	}
	
	public void init(int num, int gridSize, int iterations, double[] Is) {
		dot_width = gridSize;
		num_dots = dot_width*dot_width;
		num_vars = num;
		totalH = Is[0];
		size = 0x01 << num_vars;
		I = Is;
		D = new double[size];
		dots_per_H = (1.0-whitespace_fraction)*(double)num_dots/totalH;
		

		int[] dot_connections = new int[num_dots];
		int[][] dot_coords = new int[num_dots][2];

		
		//calculate dot_density
		for( int bit_count = num_vars; bit_count > 0; bit_count--) {
			for( int i = 0; i < size; i++) {
				if( Integer.bitCount(i) == bit_count) {
					D[i] = I[i] - totalParents(i);
				}
			}
		}

		//now attach dots
		int i2 = 1;
		double dot_remainder = 0;
		double cur_remainder = D[i2]*dots_per_H;
		for( int d = 0; d < num_dots; d++) {
			if( cur_remainder < 1) {
				dot_remainder += cur_remainder;
				if( dot_remainder >= 1) {
					dot_connections[d] = i2;
					dot_remainder--;
				} else {
					d--;
				}
				i2++;
				if( i2 >= D.length) {
					break;
				}
				cur_remainder = D[i2]*dots_per_H;
			} else {
				dot_connections[d] = i2;
				cur_remainder--;
			}
		}
		
		//now create random grids
		for( int g = 0; g < num_grids; g++) {
			int[][] grid = randomGrid(dot_connections);
			all_grids.add(new Pair<Double,int[][]>(scoreGrid(grid),grid));
		}
		Collections.sort(all_grids);
		anneal_mult = Math.exp((Math.log(2.0)-Math.log(0.01))/(double)iterations);
		for( int i = 0; i < iterations; i++) {
			int[][] grid = all_grids.get(0).b;
			perturbScored(grid,(double)dot_connections.length);
			double score = scoreGrid(grid);
			all_grids.get(0).a = score;
			if( shrink_by_area) {
				init_rate = Math.sqrt(4.0*(double)(iterations-i)/(double)iterations);
			} else {
				init_rate = 2.0*(double)(iterations-i)/(double)iterations;
			}
			//init_rate -= 2.0/(double)iterations;//*= anneal_mult;
			//init_rate *= anneal_mult;
			//System.out.println(score);
		}
		/*
		while( init_rate*(double)dot_connections.length > 1.0) {
			for( int g = 0; g < num_grids*3; g++) {
				int c = (int)(Math.random()*(double)num_grids);
				int[][] grid = cloneGrid(all_grids.get(c).b);
				perturbScored(grid,init_rate*(double)dot_connections.length);
				all_grids.add(new Pair<Double,int[][]>(scoreGrid(grid),grid));
			}
			Collections.sort(all_grids);
			while( all_grids.size() > num_grids) {
				all_grids.remove(num_grids);
			}
			init_rate *= anneal_mult;
			System.out.println(init_rate*(double)dot_connections.length+" "+all_grids.get(0).a);
		}
		*/
		
		//now print out best grid
		//System.out.println("start grid:");
		dot_grid = all_grids.get(0).b;
		for( int x = 0; x < dot_grid.length; x++) {
			for( int y = 0; y < dot_grid.length; y++) {
				int dot = dot_grid[x][y];
				String s = "";
				for( int i = 0; i < num_vars; i++) {
					s = (((0x01 << i) & dot) != 0 ? "1" : "0") + s;
				}
				System.out.print(s+",");
			}
			System.out.println();
		}
	}
	
	public int[][] randomGrid(int[] dot_connections) {
		//randomize array
		for( int d = 0; d < dot_connections.length; d++) {
			int d2 = (int)(Math.random()*(double)dot_connections.length);
			int t = dot_connections[d];
			dot_connections[d] = dot_connections[d2];
			dot_connections[d2] = t;
		}
		
		//write to grid
		int[][] dot_grid = new int[dot_width][dot_width];
		int d = 0;
		for( int x = 0; x < dot_grid.length; x++) {
			for( int y = 0; y < dot_grid.length; y++) {
				dot_grid[x][y] = dot_connections[d++];
			}
		}
		return dot_grid;
	}
	public int[][] cloneGrid(int[][] sourceGrid) {
		int[][] dot_grid = new int[dot_width][dot_width];
		for( int x = 0; x < dot_grid.length; x++) {
			for( int y = 0; y < dot_grid.length; y++) {
				dot_grid[x][y] = sourceGrid[x][y];
			}
		}
		return dot_grid;
	}
	
	public void perturbScored( int[][] dot_grid, double rate) {
		double[] counts = new double[num_vars];
		double[][] centers = new double[num_vars][2];
		
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
						double dd = dx*dx + dy*dy;
						if( !use_squared_distance) {
							dd = Math.sqrt(dd);
						}
						if( divide_by_area) {
							double area = I[(0x01 << i)];
							if( !use_squared_distance) {
								area = Math.sqrt(area);
							}
							if( area != 0) {
								dd /= area;
							}
						}
						ssd += dd; 
					}
				}
			}
		}

		int N = (int)rate; //better to poisson estimate this.
		//for( int n = 0; n < N; n++) {
		for( int x1 = 0; x1 < dot_grid.length; x1++) {
			for( int y1 = 0; y1 < dot_grid.length; y1++) {
				//int x1 = (int)((Math.random())*(double)dot_grid.length);
				//int y1 = (int)((Math.random())*(double)dot_grid.length);
				int idx = (int)((Math.random()-0.5)*init_rate*(double)dot_grid.length);
				int idy = (int)((Math.random()-0.5)*init_rate*(double)dot_grid.length);
				//int x2 = (int)((Math.random())*(double)dot_grid.length);
				//int y2 = (int)((Math.random())*(double)dot_grid.length);
				int x2 = (x1+idx) < 0 ? x1 : (x1+idx) >= dot_grid.length ? x1 : (x1+idx);
				int y2 = (y1+idy) < 0 ? y1 : (y1+idy) >= dot_grid.length ? y1 : (y1+idy);
				
				int dot1 = dot_grid[x1][y1];
				int dot2 = dot_grid[x2][y2];
				double ssd0 = 0;
				for( int i = 0; i < num_vars; i++) {
					if( ((0x01 << i) & dot1) != 0 ) {
						double dx = x1 - centers[i][0];
						double dy = y1 - centers[i][1];
						double dd = dx*dx + dy*dy;
						if( !use_squared_distance) {
							dd = Math.sqrt(dd);
						}
						if( divide_by_area) {
							double area = I[(0x01 << i)];
							if( !use_squared_distance) {
								area = Math.sqrt(area);
							}
							if( area != 0) {
								dd /= area;
							}
						}
						ssd0 -= dd;
					}
					if( ((0x01 << i) & dot2) != 0 ) {
						double dx = x2 - centers[i][0];
						double dy = y2 - centers[i][1];
						double dd = dx*dx + dy*dy;
						if( !use_squared_distance) {
							dd = Math.sqrt(dd);
						}
						if( divide_by_area) {
							double area = I[(0x01 << i)];
							if( !use_squared_distance) {
								area = Math.sqrt(area);
							}
							if( area != 0) {
								dd /= area;
							}
						}
						ssd0 -= dd;
					}
					
					if( ((0x01 << i) & dot1) != 0 ) {
						double dx = x2 - centers[i][0];
						double dy = y2 - centers[i][1];
						double dd = dx*dx + dy*dy;
						if( !use_squared_distance) {
							dd = Math.sqrt(dd);
						}
						if( divide_by_area) {
							double area = I[(0x01 << i)];
							if( !use_squared_distance) {
								area = Math.sqrt(area);
							}
							if( area != 0) {
								dd /= area;
							}
						}
						ssd0 += dd;
					}
					if( ((0x01 << i) & dot2) != 0 ) {
						double dx = x1 - centers[i][0];
						double dy = y1 - centers[i][1];
						double dd = dx*dx + dy*dy;
						if( !use_squared_distance) {
							dd = Math.sqrt(dd);
						}
						if( divide_by_area) {
							double area = I[(0x01 << i)];
							if( !use_squared_distance) {
								area = Math.sqrt(area);
							}
							if( area != 0) {
								dd /= area;
							}
						}
						ssd0 += dd;
					}
				}
				if( ssd0 < 0) {
					dot_grid[x1][y1] = dot2;
					dot_grid[x2][y2] = dot1;
				}
			}
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
						if( use_squared_distance) {
							ssd += dx*dx + dy*dy;
						} else {
							ssd += Math.sqrt(dx*dx + dy*dy);
						}
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
	public void toCsv(int[][] grid, String fn) {
		try {
			File f = new File(fn);
			FileOutputStream fis = new FileOutputStream(f);
			
			StringBuffer sb = new StringBuffer();
			for( int i = 0; i < grid.length; i++) {
				appendLine(sb,grid[i]);
			}
			fis.write(sb.toString().getBytes());
			fis.flush();
			fis.close();
			//System.out.println(sb.toString());
			System.out.println("Done.");
			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	public void appendLine(StringBuffer sb, int[] dd) {
		for( int j = 0; j < dd.length; j++) {
			int dot = dd[j];
			if( j > 0) {
				sb.append(", ");
			}
			String s = "";
			for( int i = 0; i < num_vars; i++) {
				s = (((0x01 << i) & dot) != 0 ? "1" : "0") + s;
			}
			sb.append(""+s);
		}
		sb.append("\n");
	}
	
	int[][] minusRGBs = new int[][]{
		new int[]{128,0,0},
		new int[]{0,128,0},
		new int[]{0,0,128},
		new int[]{0,64,64},
		new int[]{64,0,64},
		new int[]{64,64,0},
		
		new int[]{32,96,0},
		new int[]{0,32,96},
		new int[]{96,0,32},
		new int[]{96,32,0},
		new int[]{0,96,32},
		new int[]{32,0,96},
	};

	@Override
	public void draw(Graphics2D d, int width, int height) {
		double inc = (double)width/(double)dot_width;
		for( int x = 0; x < dot_width; x++) {
			for( int y = 0; y < dot_width; y++) {
				int x0 = (int)Math.round(inc*(double)x);
				int x1 = (int)Math.round(inc*(double)(x+1));
				int y0 = (int)Math.round(inc*(double)y);
				int y1 = (int)Math.round(inc*(double)(y+1));
				int dot = dot_grid[x][y];
				int[] c = new int[]{255,255,255};
				for( int i = 0; i < minusRGBs.length; i++) {
					if( (dot & (0x01 << i)) != 0) {
						c[0] -= minusRGBs[i][0];
						c[1] -= minusRGBs[i][1];
						c[2] -= minusRGBs[i][2];
					}
				}
				d.setColor(new Color(c[0],c[1],c[2]));
				d.fillRect(x0, y0, x1-x0, y1-y0);
			}	
		}
		
		// TODO Auto-generated method stub
		
	}

}
