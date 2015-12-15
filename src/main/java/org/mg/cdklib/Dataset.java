package org.mg.cdklib;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.mg.cdklib.Fingerprinter.Type;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.io.SMILESReader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

public class Dataset
{
	List<IAtomContainer> molecules;
	List<String> inchis;

	LinkedHashSet<String> properties;
	Set<String> uniqInchis;

	HashMap<Type, List<BitSet>> fingerprints = new HashMap<>();

	public static enum Data
	{
		Drugbank
	}

	public static Dataset get(Data data) throws FileNotFoundException, IOException, CDKException
	{
		if (data == Data.Drugbank)
			return Dataset.parseDataset("data/drugbank/all.sdf");
		throw new IllegalStateException();
	}

	public static Dataset parseDataset(String inputFile) throws FileNotFoundException, IOException, CDKException
	{
		ISimpleChemObjectReader reader;
		if (inputFile.endsWith("smi"))
			reader = new SMILESReader(new FileInputStream(inputFile));
		else
			reader = new ReaderFactory().createReader(new InputStreamReader(new FileInputStream(inputFile)));
		IChemFile content = (IChemFile) reader.read((IChemObject) new ChemFile());
		List<IAtomContainer> molecules = ChemFileManipulator.getAllAtomContainers(content);
		reader.close();
		return new Dataset(molecules);
	}

	public Dataset(Collection<IAtomContainer> molecules)
	{
		this.molecules = new ArrayList<>(molecules);
		properties = new LinkedHashSet<>();
		inchis = new ArrayList<>();
		uniqInchis = new HashSet<>();

		for (IAtomContainer mol : molecules)
			for (Object p : mol.getProperties().keySet())
				properties.add(p.toString());

		if (properties.contains("INCHI_IDENTIFIER"))
		{
			for (IAtomContainer mol : molecules)
			{
				String inchi = mol.getProperty("INCHI_IDENTIFIER").toString();
				inchis.add(inchi);
				uniqInchis.add(inchi);
			}
		}
		System.out.println("loaded dataset:\n" + this + "\n" + properties);
	}

	@Override
	public String toString()
	{
		return "#mols: " + molecules.size() + "\n#props: " + properties.size() + "\n#uniqInchis: " + uniqInchis.size();
	}

	public static final String PROB_SIM = "Similarity";

	public List<IAtomContainer> search(IAtomContainer query, Type fp, double minSim) throws CDKException
	{
		BitSet queryBitSet = Fingerprinter.get(query, fp);
		List<IAtomContainer> res = new ArrayList<>();
		int i = 0;
		for (IAtomContainer m : molecules)
		{
			double sim = Fingerprinter.tanimotoSimilarity(queryBitSet, Fingerprinter.get(m, fp));
			if (sim >= minSim)
			{
				molecules.get(i).setProperty(PROB_SIM, sim);
				res.add(molecules.get(i));
			}
			i++;
		}
		Collections.sort(res, new Comparator<IAtomContainer>()
		{
			@Override
			public int compare(IAtomContainer o1, IAtomContainer o2)
			{
				return ((Double) o2.getProperty(PROB_SIM)).compareTo((Double) o1.getProperty(PROB_SIM));
			}
		});
		return res;
	}

	public List<IAtomContainer> getMolecules()
	{
		return molecules;
	}

	//	public static void main(String[] args) throws FileNotFoundException, IOException, CDKException
	//	{
	//		Dataset drugBank = new Dataset("/home/martin/data/drugbank/all.ob.sdf");
	//		Dataset dreamCompounds = new Dataset("/home/martin/ownCloud/dream_challenge/pubchem/compounds.smi");
	//
	//		dreamCompounds.findSimilar(drugBank, 0.05, 0.2);
	//		//String sdf = "/home/martin/data/caco2/caco2.sdf";
	//		//		String sdf = "/home/martin/data/nctrer/NCTRER_v4b_232_15Feb2008.ob3d.sdf";
	//		//Dataset d = new Dataset(sdf);
	//		//		System.out.println(d.toString());
	//		//		List<IAtomContainer> sim = d.search(CDKUtil.parseSmiles("CN(CCC1)[C@@H]1C2=CC=CN=C2"), FingerprintType.ECFP6,
	//		//				0.1);
	//		//		for (IAtomContainer m : sim)
	//		//		{
	//		//			System.out.println(StringUtil.formatDouble((Double) m.getProperty(PROB_SIM)) + " " + CDKUtil.toSmiles(m));
	//		//		}
	//
	//		//		System.out.println(d.properties);
	//	}
}
