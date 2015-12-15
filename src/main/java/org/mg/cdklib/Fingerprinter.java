package org.mg.cdklib;

import java.util.BitSet;
import java.util.HashMap;

import org.mg.javalib.util.BitSetUtil;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;

public class Fingerprinter
{
	public enum Type
	{
		ECFP6;

		public int getSize()
		{
			return 1024;
		}
	}

	private static CircularFingerprinter ecfp6 = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6);
	private static HashMap<Type, HashMap<IAtomContainer, BitSet>> cache = new HashMap<>();

	public static BitSet get(IAtomContainer mol, Type fp) throws CDKException
	{
		if (!cache.containsKey(fp))
			cache.put(fp, new HashMap<IAtomContainer, BitSet>());
		HashMap<IAtomContainer, BitSet> map = cache.get(fp);
		if (!map.containsKey(mol))
		{
			BitSet bs;
			switch (fp)
			{
				case ECFP6:
					bs = ecfp6.getBitFingerprint(mol).asBitSet();
					break;
				default:
					throw new RuntimeException("not-implemented: " + fp);
			}
			map.put(mol, bs);
		}
		return map.get(mol);
	}

	public static double tanimotoSimilarity(BitSet b1, BitSet b2)
	{
		return BitSetUtil.tanimotoSimilarity(b1, b2);
	}
}
