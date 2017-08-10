
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

}
