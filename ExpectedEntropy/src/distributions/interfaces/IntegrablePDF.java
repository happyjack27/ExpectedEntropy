package distributions.interfaces;

public interface IntegrablePDF<TO,TD> {
	public double getCumulativeDistribution(TO thetas, TD value);

}
