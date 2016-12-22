package org.mg.cdklib.match;

import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.mg.cdklib.data.CDKDataset;
import org.mg.cdklib.data.DataProvider;
import org.mg.cdklib.data.DataProvider.Dataset;
import org.mg.javalib.util.StringUtil;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.KlekotaRothFingerprinter;

public class SubstructureCountTest
{
	public static void computeFPCountsForDataset() throws CDKException
	{
		SubstructureCountMatcher scm = SubstructureCountMatcherProvider
				.getInstance(new KlekotaRothFingerprinter());
		CDKDataset d = DataProvider.getDataset(Dataset.AMES);
		for (String smiles : d.getSmiles())
		{
			short[] counts = scm.matchCount(smiles);
			System.out.println(
					Arrays.stream(ArrayUtils.toObject(counts)).mapToInt(s -> s.intValue()).sum());
		}
	}

	public static void computeCount() throws CDKException
	{
		SubstructureCountMatcher scm = SubstructureCountMatcherProvider
				.getInstance(new String[] { "*" });
		//String smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";
		//		String smiles = "O1C(CO)C(O)C(O)C(O)C1O[C@@]1(OC(C(O)C1O)CO)CO";
		String smiles = "C/C=C/C[C@H]([C@H]([C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)O)C";
		short[] counts = scm.matchCount(smiles);
		for (int i = 0; i < counts.length; i++)
			if (counts[i] > 0)
				System.out.println(counts[i] + " http://localhost:8080/chem-service/depict?smiles="
						+ StringUtil.urlEncodeUTF8(smiles) + "&smarts="
						+ StringUtil.urlEncodeUTF8(scm.getSmarts(i)));
	}

	public static void main(String[] args) throws CDKException
	{
		computeCount();
	}
}
