package org.mg.cdklib.cfp;

import org.openscience.cdk.fingerprint.CircularFingerprinter;

public enum CFPType
{
	ecfp6, ecfp4, ecfp2, ecfp0, fcfp6, fcfp4, fcfp2, fcfp0;

	int getClassType()
	{
		if (this == ecfp6)
			return CircularFingerprinter.CLASS_ECFP6;
		else if (this == fcfp6)
			return CircularFingerprinter.CLASS_FCFP6;
		else if (this == ecfp4)
			return CircularFingerprinter.CLASS_ECFP4;
		else if (this == fcfp4)
			return CircularFingerprinter.CLASS_FCFP4;
		else if (this == ecfp2)
			return CircularFingerprinter.CLASS_ECFP2;
		else if (this == fcfp2)
			return CircularFingerprinter.CLASS_FCFP2;
		else if (this == ecfp0)
			return CircularFingerprinter.CLASS_ECFP0;
		else if (this == fcfp0)
			return CircularFingerprinter.CLASS_FCFP0;
		else
			throw new IllegalStateException("wtf");
	}

	public int getDiameter()
	{
		if (this == ecfp6)
			return 6;
		else if (this == fcfp6)
			return 6;
		else if (this == ecfp4)
			return 4;
		else if (this == fcfp4)
			return 4;
		else if (this == ecfp2)
			return 2;
		else if (this == fcfp2)
			return 2;
		else if (this == ecfp0)
			return 0;
		else if (this == fcfp0)
			return 0;
		else
			throw new IllegalStateException("wtf");
	}

	public boolean isECFP()
	{
		return (this == ecfp0 || this == ecfp2 || this == ecfp4 || this == ecfp6);
	}

	public String toNiceString()
	{
		return this.toString().toUpperCase();
	}
}