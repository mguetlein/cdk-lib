package org.mg.cdklib.data;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.mg.cdklib.CDKConverter;
import org.mg.javalib.util.ArrayUtil;
import org.mg.javalib.util.FileUtil;
import org.mg.javalib.util.FileUtil.CSVFile;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

public class DataLoader
{
	public static boolean VERBOSE = false;

	public static CDKDataset getDatasetFromSDF(String name, String path)
			throws CDKException, FileNotFoundException, IOException
	{
		return getDatasetFromSDF(name, path, null);
	}

	public static CDKDataset getDatasetFromSDF(String name, String path, String endpoint)
			throws CDKException, FileNotFoundException, IOException
	{
		if (VERBOSE)
			System.err.println("data-loader> " + endpoint);

		List<String> endpoints = new ArrayList<>();
		List<String> smiles = new ArrayList<>();

		ISimpleChemObjectReader reader = new ReaderFactory()
				.createReader(new InputStreamReader(new FileInputStream(path)));
		IChemFile content = (IChemFile) reader.read((IChemObject) new ChemFile());
		int invalidCompound = 0;
		int missingEndpoint = 0;
		for (IAtomContainer a : ChemFileManipulator.getAllAtomContainers(content))
		{
			if (a.getAtomCount() == 0)
				invalidCompound++;
			else if (endpoint != null && (a.getProperty(endpoint) == null
					|| a.getProperty(endpoint).toString().equals("unspecified")
					|| a.getProperty(endpoint).toString().equals("blank")
					|| a.getProperty(endpoint).toString().equals("inconclusive")))
				missingEndpoint++;
			else
			{
				String smi = new SmilesGenerator().create(a);
				//				CDKUtil.setMolForSmiles(smi, a);

				//				try
				//				{
				//					IAtomContainer m2 = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smi);
				//					if (m2.getAtomCount() == 0)
				//						throw new RuntimeException("num atoms 0");
				//					if (m2.getAtomCount() != a.getAtomCount())
				//						throw new RuntimeException("num atoms " + a.getAtomCount() + " != " + m2.getAtomCount());
				//				}
				//				catch (Exception e)
				//				{
				//					e.printStackTrace();
				//					System.exit(1);
				//				}
				smiles.add(smi);
				if (endpoint != null)
					endpoints.add(a.getProperty(endpoint).toString());
				else
					endpoints.add("n/a");
			}
		}
		reader.close();

		List<String> warnings = new ArrayList<>();
		if (invalidCompound > 0)
			warnings.add("Removed " + invalidCompound
					+ " compounds that could not be read by the CDK library.");
		if (missingEndpoint > 0)
			warnings.add("Removed " + missingEndpoint
					+ " compounds with missing/invalid enpoint values.");

		return createDataset(name, smiles, endpoints, warnings);
	}

	private static CDKDataset createDataset(String name, List<String> smiles,
			List<String> endpoints, List<String> warnings) throws CDKException
	{
		HashMap<String, HashSet<String>> uniqToActivity = new HashMap<>();
		int idx = 0;
		for (String smi : smiles)
		{
			//			System.err.println("smiles: " + smi);
			String uniq = CDKConverter.toAbsoluteSmiles(smi);
			if (!uniqToActivity.containsKey(uniq))
				uniqToActivity.put(uniq, new HashSet<String>());
			uniqToActivity.get(uniq).add(endpoints.get(idx));
			idx++;
		}

		List<String> uSmiles = new ArrayList<>();
		List<String> uEndpoints = new ArrayList<>();
		HashMap<String, Integer> uniqAdded = new HashMap<>();
		int skipEq = 0;
		int skipDiff = 0;
		for (int i = 0; i < smiles.size(); i++)
		{
			//					System.out.println(i + " " + smiles.get(i));
			String uniq = CDKConverter.toAbsoluteSmiles(smiles.get(i));
			if (uniqAdded.containsKey(uniq))
			{
				skipEq++;
				if (VERBOSE)
					System.err.println("createDataset> skip equal value " + uniqToActivity.get(uniq)
							+ " duplicate " + uniq + ":\n" + i + " " + smiles.get(i)
							+ "\nis equal to:\n" + uniqAdded.get(uniq) + " "
							+ smiles.get(uniqAdded.get(uniq)));
			}
			else if (uniqToActivity.get(uniq).size() > 1)
			{
				skipDiff++;
				if (VERBOSE)
					System.err.println("createDataset> skip different value duplicate: " + i + " "
							+ uniqToActivity.get(uniq) + " " + smiles.get(i) + " " + uniq);
			}
			else
			{
				uSmiles.add(smiles.get(i));
				uEndpoints.add(endpoints.get(i));
				uniqAdded.put(uniq, i);
			}
		}
		if (skipEq > 0)
			warnings.add("Removed " + skipEq
					+ " duplicate occurences of compounds (with equal endpoint values).");
		if (skipDiff > 0)
			warnings.add("Removed " + skipDiff
					+ " compounds that occured multiple times with different endpoint values.");

		return new CDKDataset(name, uSmiles, uEndpoints, warnings);
	}

	public static CDKDataset getDatasetFromCSV(String name, String path, int smilesCol,
			int endpointsCol) throws CDKException
	{
		List<String> endpoints = new ArrayList<>();
		List<String> smiles = new ArrayList<>();

		CSVFile csv = FileUtil.readCSV(path);
		for (int i = 1; i < csv.content.size(); i++)
		{
			smiles.add(csv.content.get(i)[smilesCol]);
			endpoints.add(csv.content.get(i)[endpointsCol]);
		}

		return createDataset(name, smiles, endpoints, new ArrayList<String>());
	}

	public static CDKDataset getDatasetFromSMILES(String name, String path) throws CDKException
	{
		return getDatasetFromSMILES(name, path, false);
	}

	public static CDKDataset getDatasetFromSMILES(String name, String path, boolean noEndpoints)
			throws CDKException
	{
		List<String> endpoints = new ArrayList<>();
		List<String> smiles = new ArrayList<>();
		int numInvalid = 0;

		for (String s : FileUtil.readStringFromFile(path).split("\n"))
		{
			String ss[] = s.split("\\s");
			try
			{
				CDKConverter.toAbsoluteSmiles(ss[0]);
				smiles.add(ss[0]);
				endpoints.add((ss.length > 1 && !noEndpoints) ? ss[1] : "");
			}
			catch (CDKException e)
			{
				numInvalid++;
				if (VERBOSE)
					System.err.println("createDataset> could not read smiles: " + ss[0]);
			}
		}

		List<String> warnings = new ArrayList<>();
		if (numInvalid > 0)
			warnings.add("Could not read " + numInvalid + " smiles.");
		return createDataset(name, smiles, endpoints, new ArrayList<String>());
	}

	public static String[] getClassValues(List<String> endpoints)
	{
		String classValues[] = ArrayUtil.removeDuplicates(ArrayUtil.toArray(endpoints));
		Arrays.sort(classValues);
		return classValues;
	}

	public static Integer getActiveIdx(String[] classValues)
	{
		Integer activeIdx = null;
		for (int i = 0; i < classValues.length; i++)
			if (classValues[i].equals("active") || classValues[i].equals("mutagen")
					|| classValues[i].equals("1") || classValues[i].equals("most-concern"))
				activeIdx = i;
		if (activeIdx == null)
			throw new IllegalStateException("what is active? " + ArrayUtil.toString(classValues));
		return activeIdx;
	}

}
