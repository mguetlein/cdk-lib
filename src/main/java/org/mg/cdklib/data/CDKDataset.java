package org.mg.cdklib.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.mg.cdklib.CDKConverter;
import org.openscience.cdk.exception.CDKException;

public class CDKDataset
{
	protected List<String> smiles;
	protected List<String> endpoints;
	protected List<String> warnings;

	public CDKDataset(List<String> smiles, List<String> endpoints, List<String> warnings)
	{

		try
		{
			HashMap<String, HashSet<String>> inchiToActivity = new HashMap<>();
			int idx = 0;
			for (String smi : smiles)
			{
				String inch = CDKConverter.toInchi(smi);
				if (!inchiToActivity.containsKey(inch))
					inchiToActivity.put(inch, new HashSet<String>());
				inchiToActivity.get(inch).add(endpoints.get(idx));
				idx++;
			}

			List<String> uSmiles = new ArrayList<>();
			List<String> uEndpoints = new ArrayList<>();
			HashSet<String> inchiAdded = new HashSet<>();
			int skipEq = 0;
			int skipDiff = 0;
			for (int i = 0; i < smiles.size(); i++)
			{
				//					System.out.println(i + " " + smiles.get(i));
				String inchi = CDKConverter.toInchi(smiles.get(i));
				if (inchiAdded.contains(inchi))
				{
					skipEq++;
					//						System.out.println("skip equal value duplicate: " + i + " " + inchiToActivity.get(inchi) + " "
					//								+ smiles.get(i) + " " + inchi);
				}
				else if (inchiToActivity.get(inchi).size() > 1)
				{
					skipDiff++;
					//						System.out.println("skip different value duplicate: " + i + " " + inchiToActivity.get(inchi)
					//								+ " " + smiles.get(i) + " " + inchi);
				}
				else
				{
					uSmiles.add(smiles.get(i));
					uEndpoints.add(endpoints.get(i));
					inchiAdded.add(inchi);
				}
			}
			if (skipEq > 0)
				warnings.add("Removed " + skipEq + " duplicate occurences of compounds (with equal endpoint values).");
			if (skipDiff > 0)
				warnings.add("Removed " + skipDiff
						+ " compounds that occured multiple times with different endpoint values.");

			this.smiles = uSmiles;
			this.endpoints = uEndpoints;
			this.warnings = warnings;
		}
		catch (CDKException e)
		{
			throw new RuntimeException(e);
		}
	}

	public List<String> getEndpoints()
	{
		return endpoints;
	}

	public List<String> getSmiles()
	{
		return smiles;
	}

	public List<String> getWarnings()
	{
		return warnings;
	}
}