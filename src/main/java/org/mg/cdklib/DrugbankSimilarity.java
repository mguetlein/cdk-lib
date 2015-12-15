package org.mg.cdklib;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;

import org.mg.cdklib.Fingerprinter.Type;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

public class DrugbankSimilarity
{
	public static HashMap<String, HashMap<String, Double>> get(Dataset data, String propKey, double relMinFreq,
			double simThres) throws FileNotFoundException, IOException, CDKException
	{
		int absMinFreq = Math.max(1, (int) (data.getMolecules().size() * relMinFreq));
		System.err.println("looking for drug-bank compounds with sim>=" + simThres + " in at least " + absMinFreq
				+ " compounds");

		Dataset drugbank = Dataset.get(Dataset.Data.Drugbank);
		List<IAtomContainer> selectedDrugbankCompounds = new ArrayList<>();
		for (IAtomContainer m : drugbank.getMolecules())
		{
			List<IAtomContainer> similarCompounds = data.search(m, Type.ECFP6, simThres);
			if (similarCompounds.size() >= absMinFreq)
			{
				//				IAtomContainer[] draw = new IAtomContainer[1 + similarCompounds.size()];
				//				m.setProperty("Sim", "n/a");
				//				draw[0] = m;
				//				for (int i = 0; i < similarCompounds.size(); i++)
				//				{
				//					similarCompounds.get(i).setProperty("Sim",
				//							StringUtil.formatDouble((Double) similarCompounds.get(i).getProperty(PROB_SIM)));
				//					draw[i + 1] = similarCompounds.get(i);
				//				}
				//				Graphics.draw(draw, ListUtil.createList("Sim", "GENERIC_NAME"));
				//				System.out.println();

				selectedDrugbankCompounds.add(m);
			}
		}

		System.err.println("found " + selectedDrugbankCompounds.size() + " drug-bank compounds");

		HashMap<String, HashMap<String, Double>> res = new HashMap<>();
		for (IAtomContainer m : selectedDrugbankCompounds)
		{
			String name = m.getProperty("GENERIC_NAME");
			HashMap<String, Double> sim = new HashMap<>();
			BitSet mBS = Fingerprinter.get(m, Fingerprinter.Type.ECFP6);
			for (IAtomContainer mol : data.getMolecules())
			{
				String k = mol.getProperty(propKey);
				double d = Fingerprinter.tanimotoSimilarity(mBS, Fingerprinter.get(mol, Fingerprinter.Type.ECFP6));
				sim.put(k, d);
			}
			res.put(name, sim);
		}
		return res;
	}
}
