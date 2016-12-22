package org.mg.cdklib.match;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.mg.cdklib.CDKConverter;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

public class DefaultSubstructureCountMatcher implements SubstructureCountMatcher, Serializable
{
	private static final long VERSION = 9L;
	private static final long serialVersionUID = VERSION;

	transient List<SMARTSQueryTool> queryTools;
	List<String> smarts;
	Map<String, BitSet> smilesToMatch;
	Map<String, short[]> smilesToCount;

	DefaultSubstructureCountMatcher(String smarts[])
	{
		this.smarts = new ArrayList<>();
		for (String smrt : smarts)
			this.smarts.add(smrt);
		smilesToMatch = new HashMap<>();
		smilesToCount = new HashMap<>();
	}

	public static long getVersion()
	{
		return VERSION;
	}

	private void doMatch(String smiles) throws CDKException
	{
		if (queryTools == null)
			queryTools = new ArrayList<>();
		IAtomContainer mol = CDKConverter.parseSmiles(smiles);
		BitSet resMatch = new BitSet(smarts.size());
		short[] resCount = new short[smarts.size()];
		for (int j = 0; j < smarts.size(); j++)
		{
			if (queryTools.size() < j + 1)
				queryTools.add(
						new SMARTSQueryTool(smarts.get(j), SilentChemObjectBuilder.getInstance()));
			if (queryTools.get(j).matches(mol))
			{
				resMatch.set(j);
				resCount[j] = (short) queryTools.get(j).getUniqueMatchingAtoms().size();
			}
		}
		smilesToMatch.put(smiles, resMatch);
		smilesToCount.put(smiles, resCount);
	}

	public short[] matchCount(String smiles) throws CDKException
	{
		if (!smilesToCount.containsKey(smiles))
			doMatch(smiles);
		return smilesToCount.get(smiles);
	}

	public BitSet match(String smiles) throws CDKException
	{
		if (!smilesToMatch.containsKey(smiles))
			doMatch(smiles);
		return smilesToMatch.get(smiles);
	}

	public int numMatchedSmiles()
	{
		return smilesToCount.size();
	}

	@Override
	public int size()
	{
		return smarts.size();
	}

	@Override
	public String getSmarts(int idx)
	{
		return smarts.get(idx);
	}

}
