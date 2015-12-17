package org.mg.cdklib.data;

import java.io.Serializable;
import java.util.List;

import org.mg.javalib.util.HashUtil;

public class CDKDataset implements Serializable
{
	private static final long serialVersionUID = 2L;

	protected String datasetName;
	protected List<String> smiles;
	protected List<String> endpoints;
	protected List<String> warnings;

	public CDKDataset(String datasetName, List<String> smiles, List<String> endpoints, List<String> warnings)
	{
		this.smiles = smiles;
		this.datasetName = datasetName;
		this.endpoints = endpoints;
		this.warnings = warnings;
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

	public String getDatasetName()
	{
		return datasetName;
	}

	@Override
	public int hashCode()
	{
		return HashUtil.hashCode(smiles, datasetName, endpoints, warnings);
	}

	@Override
	public boolean equals(Object obj)
	{
		return (obj instanceof CDKDataset //
				&& ((CDKDataset) obj).getSmiles().equals(smiles) //
				&& ((CDKDataset) obj).getDatasetName().equals(datasetName) //
				&& ((CDKDataset) obj).getEndpoints().equals(endpoints) //
		&& ((CDKDataset) obj).getWarnings().equals(warnings));
	}
}