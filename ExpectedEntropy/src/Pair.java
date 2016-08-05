
public class Pair<A extends Comparable<A>,B> implements Comparable< Pair<A,B>> {
	public A a;
	public B b;
	
	public Pair(A _a, B _b) {
		super();
		a = _a;
		b = _b;
	}
	@Override
	public int compareTo(Pair<A, B> o) {
		return a.compareTo(o.a);
	}

}
