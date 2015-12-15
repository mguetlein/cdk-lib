package org.mg.cdklib.cfp;

public enum FeatureSelection
{
	filt, fold, none;

	public String toNiceString()
	{
		switch (this)
		{
			case filt:
				return "Filtering";
			case fold:
				return "Folding";
			case none:
				return "Unprocessed";
		}
		throw new IllegalStateException();
	}

	public String attribute()
	{
		switch (this)
		{
			case filt:
				return "Filtered";
			case fold:
				return "Folded";
			case none:
				return "Unprocessed";
		}
		throw new IllegalStateException();
	}

	public String toNiceShortString()
	{
		switch (this)
		{
			case filt:
				return "Filt.";
			case fold:
				return "Fold.";
			case none:
				return "Unproc.";
		}
		throw new IllegalStateException();
	}
}