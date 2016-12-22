package org.mg.cdklib.match;

import java.util.BitSet;

import org.openscience.cdk.exception.CDKException;

public interface SubstructureCountMatcher
{
	short[] matchCount(String smiles) throws CDKException;

	BitSet match(String smiles) throws CDKException;

	int size();

	String getSmarts(int idx);
}
