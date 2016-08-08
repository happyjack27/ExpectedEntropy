package distributions;

public interface hasMarginalLikelihood<TSupport, TData> {
	public TSupport getMarginalLikelihood(TData data); 

}
