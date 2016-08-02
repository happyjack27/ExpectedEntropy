package distributions.interfaces;

public interface hasMLE<TEstimator,TData> {
	public TEstimator getMLE(TData[] data);
}
